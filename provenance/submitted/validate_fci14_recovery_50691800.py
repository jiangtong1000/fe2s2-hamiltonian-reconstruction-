"""Compare independent CAS14 FCI with DMRG and its full-space embedding.

Run in a compute allocation after both states and direct molecular checks finish.
This does not replace the external-virtual validation in validate_full.py.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import numpy as np
from pyscf import lib

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'build'))
from integrals import embedded_singlet_energy

def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    path=ROOT/'cas14';fci_path=path/'fci_validation'
    fm=json.loads((fci_path/'energies.json').read_text())
    dm=json.loads((path/'energies.json').read_text())
    digest=hashlib.sha256((path/'FCIDUMP').read_bytes()).hexdigest()
    assert fm['converged'] and fm['fcidump_sha256']==dm['fcidump_sha256']==digest
    stages=json.loads((path/'convergence.json').read_text())
    assert abs(stages[-1]['energy_total']-dm['energy_total'])<1e-7
    d1=np.load(fci_path/'rdm1.npy');d2=np.load(fci_path/'rdm2_block2_order.npy')
    h=np.load(path/'h1.npy');g=np.load(path/'eri.npy');e0=float((path/'ecore.txt').read_text())
    cas_energy=float(e0+np.einsum('pq,pq',h,d1)+.5*np.einsum('pqrs,prsq',g,d2))
    maps=json.loads((ROOT/'full_hamiltonian/index_maps.json').read_text())
    embedded=embedded_singlet_energy(np.load(ROOT/'full_hamiltonian/h1.npy'),
        np.load(ROOT/'full_hamiltonian/chol.npy',mmap_mode='r'),
        float((ROOT/'full_hamiltonian/h0.txt').read_text()),
        maps['I']+maps['K6'],maps['A14'],d1,d2)
    direct=np.load(ROOT/'validation/cas14_direct_integrals.npz')
    direct_energy=float(direct['ecore']+np.einsum('pq,pq',direct['h1'],d1)+
                        .5*np.einsum('pqrs,prsq',direct['eri'],d2))
    delta=float(dm['energy_total']-fm['energy_total'])
    occupations_fci=np.linalg.eigvalsh(d1)[::-1]
    occupations_dmrg=np.linalg.eigvalsh(np.load(path/'rdm1.npy'))[::-1]
    checks={
        'fci_rdm_energy':{'error_Eh':cas_energy-fm['energy_total'],
                          'passed':abs(cas_energy-fm['energy_total'])<1e-7},
        'fci_full_embedding':{'error_Eh':embedded-fm['energy_total'],
                              'passed':abs(embedded-fm['energy_total'])<1e-7},
        'fci_direct_AO_energy':{'error_Eh':direct_energy-fm['energy_total'],
                               'passed':abs(direct_energy-fm['energy_total'])<1e-5},
        'dmrg_vs_fci':{'dmrg_minus_fci_Eh':delta,'passed':-1e-7<=delta<1e-5},
        'rdm_traces':{'one_body':float(np.trace(d1)),
                      'two_body':float(np.einsum('ijji',d2)),
                      'passed':bool(abs(np.trace(d1)-18)<1e-7 and abs(np.einsum('ijji',d2)-306)<1e-5)}}
    report=dict(status='passed' if all(x['passed'] for x in checks.values()) else 'failed',
        checks=checks,fcidump_sha256=digest,fci_energy_total=fm['energy_total'],
        dmrg_energy_total=dm['energy_total'],fci_residual_norm=fm['residual_norm'],
        sorted_occupations_fci=occupations_fci.tolist(),sorted_occupations_dmrg=occupations_dmrg.tolist(),
        maximum_occupation_difference=float(abs(occupations_fci-occupations_dmrg).max()),
        fci_generating_job=fm['generating_slurm_job_id'],dmrg_generating_job=dm['generating_slurm_job_id'],
        validation_generating_job=os.environ.get('SLURM_JOB_ID'),
        interpretation='Independent iterative FCI ground-state check; residual and solver convergence remain part of the evidence.')
    (ROOT/'validation/cas14_fci_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    assert report['status']=='passed', 'CAS14 FCI validation failed; inspect recorded metrics'

if __name__=='__main__':main()
