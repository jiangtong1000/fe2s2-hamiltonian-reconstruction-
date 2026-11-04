"""Refresh state checks while retaining independently validated fixed integrals."""
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
from cross_basis_energy import evaluate as cross_basis_energy

def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):digest.update(block)
    return digest.hexdigest()

def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    baseline=ROOT/'validation/initial_production_50707994'
    report=json.loads((baseline/'validation_metrics.json').read_text())
    assert report['status']=='failed'
    expected={'cas20_variational_convergence_crosscheck','cas20_solver_convergence'}
    assert set(report['failed_checks'])==expected, 'Unresolved non-solver baseline failures'
    hashes=json.loads((baseline/'input_hashes.json').read_text())['inputs']
    # Fresh streaming hashes include the complete 26-GiB factor file. A changed
    # Hamiltonian or orbital basis invalidates reuse of the molecular checks.
    for name,record in hashes.items():
        path=ROOT/name
        assert path.stat().st_size==record['size_bytes'],name
        assert sha(path)==record['sha256'],f'Validated input changed: {name}'
    checks=report['checks']
    cross=cross_basis_energy()
    checks['cas20_cross_basis_trial']=cross
    checks['cas20_variational_convergence_crosscheck']=dict(
        energy_error_lower_bound=cross['production_energy_error_lower_bound'],
        passed=cross['production_within_1e_5_of_saved_trial'],
        note='Necessary check against an existing trial; not sufficient to certify the ground state')
    h=np.load(ROOT/'full_hamiltonian/h1.npy')
    l=np.load(ROOT/'full_hamiltonian/chol.npy',mmap_mode='r')
    e0=float((ROOT/'full_hamiltonian/h0.txt').read_text())
    maps=json.loads((ROOT/'full_hamiltonian/index_maps.json').read_text())
    states={}
    for name,core,active,ne in [('cas20',maps['I'],maps['cas20_local_to_full_zero_based'],30),
                               ('cas14',maps['I']+maps['K6'],maps['A14'],18)]:
        path=ROOT/name
        meta=json.loads((path/'energies.json').read_text())
        stages=json.loads((path/'convergence.json').read_text())
        d1=np.load(path/'rdm1.npy');d2=np.load(path/'rdm2.npy')
        assert meta['fcidump_sha256']==sha(path/'FCIDUMP')
        assert abs(meta['energy_total']-stages[-1]['energy_total'])<1e-7
        assert abs(d1-np.load(path/f"rdm1_stage{stages[-1]['stage']}.npy")).max()<1e-8
        assert abs(np.trace(d1)-ne)<1e-7 and abs(np.einsum('ijji',d2)-ne*(ne-1))<1e-5
        assert abs(meta['norm']-1)<1e-8 and abs(meta['spin_square'])<1e-8
        assert abs(meta['native_state_reload_validation']['energy_error'])<1e-7
        assert abs(meta['native_state_reload_validation']['norm']-1)<1e-8
        energy=meta['energy_total']
        embedded=embedded_singlet_energy(h,l,e0,core,active,d1,d2)
        direct=np.load(baseline/f'{name}_direct_integrals.npz')
        def expectation(hh,gg,ee):
            return float(ee+np.einsum('pq,pq',hh,d1)+.5*np.einsum('pqrs,prsq',gg,d2))
        direct_energy=expectation(direct['h1'],direct['eri'],direct['ecore'])
        checks[name+'_embedded_and_direct_energies']=dict(energy_solver=energy,
            energy_embedded=embedded,energy_direct=direct_energy,embedded_error=embedded-energy,
            direct_error=direct_energy-energy,passed=abs(embedded-energy)<1e-7 and abs(direct_energy-energy)<1e-5)
        change=float(stages[-1]['energy_total']-stages[-2]['energy_total'])
        checks[name+'_solver_convergence']=dict(last_bond_stage_energy_change=change,
            target_Eh=1e-5,passed=abs(change)<1e-5,
            caveat='Successive stages provide convergence evidence, not a rigorous absolute-error bound')
        tight=expectation(np.load(path/'h1.npy'),np.load(path/'eri.npy'),float((path/'ecore.txt').read_text()))
        loose_path=ROOT/(name+'_loose')
        loose=expectation(np.load(loose_path/'h1.npy'),np.load(loose_path/'eri.npy'),float((loose_path/'ecore.txt').read_text()))
        assert abs(tight-energy)<1e-7
        checks[name+'_factorization_sensitivity']=dict(tight_state_tight_energy=tight,
            tight_state_loose_energy=loose,loose_minus_tight=loose-tight,passed=abs(loose-tight)<1e-5,
            note='Fixed-state expectation sensitivity, not a bound on solver error')
        states[name]=dict(generating_job=meta['generating_slurm_job_id'],energy_total=energy,
            hashes={n:sha(path/n) for n in ['energies.json','convergence.json','rdm1.npy','rdm2.npy']})
    failed=[k for k,v in checks.items() if v.get('passed') is False]
    report.update(status='failed' if failed else 'numerical_checks_passed',failed_checks=failed,
        fixed_integral_checks_reused_from='validation/initial_production_50707994/validation_metrics.json',
        baseline_report_sha256=sha(baseline/'validation_metrics.json'),
        validated_inputs_manifest_sha256=sha(baseline/'input_hashes.json'),
        reused_inputs_verified_by_fresh_SHA256=True,state_provenance=states,
        state_validation_job=os.environ.get('SLURM_JOB_ID'))
    (ROOT/'validation/validation_metrics.json').write_text(json.dumps(report,indent=2)+'\n')
    print('REFINED_STATE_VALIDATION',report['status'],failed,flush=True)
    assert not failed, 'Refined state checks failed; inspect validation_metrics.json'

if __name__=='__main__':main()
