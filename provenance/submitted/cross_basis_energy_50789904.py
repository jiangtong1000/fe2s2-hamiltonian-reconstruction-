"""Evaluate the saved provisional CAS20 state in the production Hamiltonian."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def evaluate():
    provisional=ROOT/'cas20_provisional'
    meta=json.loads((provisional/'energies.json').read_text())
    u=np.load(ROOT/'orbitals/U_natural_20.npy')
    assert abs(u.T@u-np.eye(20)).max()<1e-10
    ci=np.load(ROOT/'orbitals/mo_coeff_20_initial.npy')
    cf=np.load(ROOT/'orbitals/mo_coeff_20_final.npy')
    assert abs(ci@u-cf).max()<1e-10
    d1=np.load(provisional/'rdm1.npy');d2=np.load(provisional/'rdm2.npy')
    h=np.load(ROOT/'cas20/h1.npy');g=np.load(ROOT/'cas20/eri.npy')
    e0=float((ROOT/'cas20/ecore.txt').read_text())
    hi=u@h@u.T
    gi=np.einsum('pa,qb,rc,sd,abcd->pqrs',u,u,u,u,g,optimize=True)
    energy=float(e0+np.einsum('pq,pq',hi,d1)+.5*np.einsum('pqrs,prsq',gi,d2))
    # Independently rotate densities and compare the same expectation.
    df=u.T@d1@u
    gf=np.einsum('pa,rb,sc,qd,prsq->abcd',u,u,u,u,d2,optimize=True)
    energy_rotated=float(e0+np.einsum('pq,pq',h,df)+.5*np.einsum('pqrs,prsq',g,gf))
    assert abs(energy-energy_rotated)<1e-9
    result=dict(provisional_state_production_energy=energy,
        provisional_state_direct_energy=meta['energy_total'],
        production_minus_direct_same_state=energy-meta['energy_total'],
        integral_vs_density_rotation_error=energy-energy_rotated,
        production_fcidump_sha256=hashlib.sha256((ROOT/'cas20/FCIDUMP').read_bytes()).hexdigest(),
        source_rdm_sha256={n:hashlib.sha256((provisional/n).read_bytes()).hexdigest() for n in ['rdm1.npy','rdm2.npy']},
        source_state_job=meta['generating_slurm_job_id'],
        interpretation='Variational trial energy in the same production Hamiltonian. A higher production DMRG energy implies at least that energy error; this trial is not an exact ground-state bound from below.',
        passed=abs(energy-meta['energy_total'])<1e-5)
    final=ROOT/'cas20/energies.json'
    if final.exists():
        fm=json.loads(final.read_text())
        assert fm['fcidump_sha256']==result['production_fcidump_sha256']
        excess=float(fm['energy_total']-energy)
        result.update(production_energy_excess_above_saved_trial=excess,
            production_energy_error_lower_bound=max(0.,excess),
            production_within_1e_5_of_saved_trial=bool(excess<=1e-5),
            production_state_job=fm['generating_slurm_job_id'])
    (ROOT/'validation/cas20_cross_basis_energy.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(evaluate(),indent=2))
