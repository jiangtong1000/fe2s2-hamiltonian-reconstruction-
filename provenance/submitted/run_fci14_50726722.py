"""Independent singlet CAS14 FCI check of the production Hamiltonian."""
from pathlib import Path
import hashlib
import json
import math
import os
import shutil
import time
import numpy as np
from pyscf import fci,lib

ROOT=Path(__file__).resolve().parents[1]
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    path=ROOT/'cas14';out=path/'fci_validation';out.mkdir(exist_ok=True)
    h=np.load(path/'h1.npy');g=np.load(path/'eri.npy');e0=float((path/'ecore.txt').read_text())
    assert h.shape==(14,14) and g.shape==(14,)*4
    solver=fci.direct_spin0.FCI()
    solver.conv_tol=1e-10;solver.max_cycle=400;solver.max_space=64
    solver.max_memory=24000;solver.verbose=5
    ci0=None
    if (out/'ci.npy').exists():
        ci0=np.load(out/'ci.npy')
        shutil.copy2(out/'ci.npy',out/f"ci_before_job_{os.environ.get('SLURM_JOB_ID','local')}.npy")
    start=time.time()
    energy,ci=solver.kernel(h,g,14,(9,9),ci0=ci0,ecore=e0)
    np.save(out/'ci.npy',np.asarray(ci))
    h2=solver.absorb_h1e(h,g,14,(9,9),.5)
    hc=solver.contract_2e(h2,ci,14,(9,9))
    residual=float(np.linalg.norm(np.asarray(hc)-(energy-e0)*np.asarray(ci)))
    (out/'checkpoint.json').write_text(json.dumps(dict(energy_total=float(energy),converged=bool(solver.converged),
        residual_norm=residual,generating_slurm_job_id=os.environ.get('SLURM_JOB_ID')),indent=2)+'\n')
    assert solver.converged,'FCI Davidson did not converge'
    dm1,dm2=solver.make_rdm12(ci,14,(9,9))
    rdm_energy=float(e0+np.einsum('pq,pq',h,dm1)+.5*np.einsum('pqrs,pqrs',g,dm2))
    ss,multiplicity=fci.spin_op.spin_square(ci,14,(9,9))
    norm=float(np.vdot(ci,ci).real)
    assert abs(norm-1)<1e-8 and abs(ss)<1e-7,(norm,ss)
    assert abs(rdm_energy-energy)<1e-7,(rdm_energy,energy)
    np.save(out/'rdm1.npy',dm1)
    np.save(out/'rdm2_block2_order.npy',dm2.transpose(0,2,3,1))
    strings=fci.cistring.make_strings(range(14),9)
    np.save(out/'determinant_strings.npy',strings)
    data=dict(solver='PySCF direct_spin0 FCI',converged=bool(solver.converged),energy_total=float(energy),
        energy_active=float(energy-e0),ecore=e0,energy_rdm=rdm_energy,norm=norm,spin_square=float(ss),
        norb=14,nalpha=9,nbeta=9,determinant_dimension=math.comb(14,9)**2,
        convergence_tolerance=solver.conv_tol,elapsed_seconds=time.time()-start,
        residual_norm=residual,
        generating_slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        fcidump_sha256=hashlib.sha256((path/'FCIDUMP').read_bytes()).hexdigest(),
        coefficient_order='ci[alpha_string_index,beta_string_index]; ascending binary bit strings in determinant_strings.npy; bit i occupies spatial orbital i',
        spin_orbital_order='all local alpha orbitals ascending, then all local beta orbitals ascending',
        rdm2_convention='rdm2_block2_order[p,r,s,q] = sum_st <p_s^dagger r_t^dagger s_t q_s>')
    (out/'energies.json').write_text(json.dumps(data,indent=2)+'\n')
    print('CAS14_FCI_COMPLETE',json.dumps(data),flush=True)

if __name__=='__main__':main()
