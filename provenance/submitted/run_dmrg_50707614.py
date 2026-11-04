"""Identical singlet DMRG/RDM analysis for reference and reconstructed FCIDUMPs.

Run with the directory containing FCIDUMP as cwd. Retain MPS in native_state/.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import os
import numpy as np
from pyscf import ao2mo
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rdm-only',action='store_true',help='Reuse saved GS MPS without optimization')
    args=ap.parse_args()
    out=Path.cwd()
    nthreads=int(os.environ.get('OMP_NUM_THREADS','1'))
    driver=DMRGDriver(scratch=str(out/'native_state'),symm_type=SymmetryTypes.SU2,
                      n_threads=nthreads,mpi=False,stack_mem=int(os.environ.get('DMRG_STACK_GB','32'))*1024**3)
    driver.bw.b.Random.rand_seed(10405)
    driver.read_fcidump(str(out/'FCIDUMP'))
    n,ne=driver.n_sites,driver.n_elec
    driver.initialize_system(n_sites=n,n_elec=ne,spin=0,orb_sym=None)
    mpo=driver.get_qc_mpo(h1e=driver.h1e,g2e=driver.g2e,ecore=driver.ecore,
                          algo_type=MPOAlgorithmTypes.Conventional,iprint=1)
    occ_file=out/'initial_occupations.npy'
    occs=np.load(occ_file).tolist() if occ_file.exists() else None
    ket=driver.load_mps('GS') if args.rdm_only else driver.get_random_mps(tag='GS',bond_dim=300,nroots=1,occs=occs)
    results=json.loads((out/'convergence.json').read_text()) if args.rdm_only else []
    seed_path=out/'optimization_provenance.json'
    if not args.rdm_only:
        seed_path.write_text(json.dumps({'random_seed':10405},indent=2)+'\n')
    optimization_provenance=json.loads(seed_path.read_text()) if seed_path.exists() else {'random_seed':None,'note':'Original optimization seed was not recorded; RDM recovery seed does not establish it.'}
    prev_occ=None
    stages=[] if args.rdm_only else [[300]*2+[500]*2+[800]*2+[1200]*6,[1600]*6]
    for stage, dims in enumerate(stages):
        noises=([1e-4]*2+[1e-5]*2+[1e-6]*2+[0.]*6) if stage==0 else [0.]*6
        energy=driver.dmrg(mpo,ket,n_sweeps=len(dims),bond_dims=dims,noises=noises,
                           thrds=[1e-9]*len(dims),tol=1e-8,iprint=2,
                           twosite_to_onesite=len(dims)-2)
        rdm1=driver.get_1pdm(ket,iprint=1)
        occ=np.linalg.eigvalsh(rdm1)[::-1]
        assert abs(np.trace(rdm1)-ne)<1e-7
        assert np.max(abs(rdm1-rdm1.T))<1e-8
        np.save(out/f'rdm1_stage{stage}.npy',rdm1)
        result=dict(stage=stage,max_bond_dimension=max(dims),energy_total=float(energy),
                    occupations=occ.tolist(),occupation_max_change=None if prev_occ is None else float(abs(occ-prev_occ).max()))
        results.append(result)
        (out/'convergence.json').write_text(json.dumps(results,indent=2)+'\n')
        prev_occ=occ
    if args.rdm_only:
        energy=results[-1]['energy_total']
        rdm1=driver.get_1pdm(ket,iprint=1)
        occ=np.linalg.eigvalsh(rdm1)[::-1]
    rdm2=driver.get_2pdm(ket,iprint=1)
    # block2 dm2[p,r,s,q] = sum_st <p_s^+ r_t^+ s_t q_s>.
    eri4=ao2mo.restore(1,driver.g2e,n)
    np.save(out/'rdm1.npy',rdm1)
    np.save(out/'rdm2.npy',rdm2)
    energy_rdm=float(driver.ecore+np.einsum('pq,pq',driver.h1e,rdm1)+.5*np.einsum('pqrs,prsq',eri4,rdm2))
    assert abs(energy_rdm-energy)<1e-7,(energy_rdm,energy)
    assert abs(np.einsum('ijji',rdm2)-ne*(ne-1))<1e-5
    np.save(out/'rdm1.npy',rdm1)
    np.save(out/'rdm2.npy',rdm2)
    np.savetxt(out/'occupations.txt',occ,fmt='%.14f')
    identity=driver.get_identity_mpo()
    norm=float(driver.expectation(ket,identity,ket))
    spin_square=float(driver.expectation(ket,driver.get_spin_square_mpo(),ket))/norm
    assert abs(norm-1)<1e-8,norm
    assert abs(spin_square)<1e-8,spin_square
    meta=dict(norb=n,nelec=ne,spin=0,symmetry='SU2 singlet',energy_total=float(energy),
              ecore=float(driver.ecore),energy_active=float(energy-driver.ecore),energy_rdm=energy_rdm,norm=norm,
              fcidump_sha256='',spin_square=spin_square,random_seed=optimization_provenance['random_seed'])
    meta['optimization_provenance']=optimization_provenance
    meta['rdm_only_recovery']=args.rdm_only
    meta['fcidump_sha256']=hashlib.sha256((out/'FCIDUMP').read_bytes()).hexdigest()
    meta['block2_version']=importlib.metadata.version('block2-mpi')
    meta['solver_site_to_local_zero_based']=list(range(n))
    meta['rdm1_convention']='spin-summed <a_p^dagger a_q>'
    meta['rdm2_convention']='dm2[p,r,s,q] = sum_sigma,tau <a_p,sigma^dagger a_r,tau^dagger a_s,tau a_q,sigma>'
    meta['native_state']=dict(directory='native_state',tag='GS',spin_adapted=True)
    (out/'energies.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('DMRG_RDM_COMPLETE',json.dumps(meta))
if __name__=='__main__': main()
