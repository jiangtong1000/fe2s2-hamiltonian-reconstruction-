"""Exact MPS expectation of a product of local double-occupation projectors.

For SU2, the same scalar four-fermion expression used by block2 for g[i,i,i,i]
with coefficient one is n_i,alpha n_i,beta. Scalar products couple with unit CG.
"""
from pathlib import Path
import argparse
import json
import os
import hashlib
import numpy as np
from pyblock2.driver.core import DMRGDriver,SymmetryTypes,MPOAlgorithmTypes

ROOT=Path(__file__).resolve().parents[1]

def projector_mpo(driver,sites):
    assert sites==sorted(set(sites)) and sites
    local='((C+(C+D)0)1+D)0'
    expr=local
    for _ in sites[1:]:expr=f'({expr}+{local})0'
    builder=driver.expr_builder()
    builder.add_term(expr,[i for i in sites for _ in range(4)],1.0)
    # Indices are already sorted. Automatic SU2 operator-order expansion at
    # length 24 is unnecessary and can exhaust memory.
    return driver.get_mpo(builder.finalize(adjust_order=False),iprint=0)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    path=ROOT/'validation/sector_projector_smoke' if args.smoke else ROOT/'cas20'
    path.mkdir(exist_ok=True)
    if args.smoke:
        from pyscf import gto,scf,ao2mo,tools
        mol=gto.M(atom='H 0 0 0; H 0 0 1.1; H 0 0 2.5; H 0 0 4.2',basis='sto-3g',verbose=0)
        mf=scf.RHF(mol).run()
        h=mf.mo_coeff.T@mf.get_hcore()@mf.mo_coeff;g=ao2mo.kernel(mol,mf.mo_coeff)
        tools.fcidump.from_integrals(str(path/'FCIDUMP'),h,g,4,4,nuc=mol.energy_nuc())
    # Separate contraction scratch preserves the source MPS files.
    import shutil
    scratch=path/'sector_scratch';scratch.mkdir(exist_ok=True)
    for p in ([] if args.smoke else (path/'native_state').glob('*')):
        if p.is_file() and ('MPS.GS.' in p.name or 'MPS.INFO.GS.' in p.name or p.name=='GS-mps_info.bin'):
            shutil.copy2(p,scratch/p.name)
    driver=DMRGDriver(scratch=str(scratch),symm_type=SymmetryTypes.SU2,
        n_threads=int(os.environ.get('OMP_NUM_THREADS','1')),stack_mem=4*1024**3,mpi=False)
    driver.read_fcidump(str(path/'FCIDUMP'))
    driver.initialize_system(n_sites=driver.n_sites,n_elec=driver.n_elec,spin=0,orb_sym=None)
    if args.smoke:
        ket=driver.get_random_mps(tag='GS',bond_dim=32,nroots=1)
        mpo=driver.get_qc_mpo(h1e=driver.h1e,g2e=driver.g2e,ecore=driver.ecore,algo_type=MPOAlgorithmTypes.Conventional,iprint=0)
        driver.dmrg(mpo,ket,n_sweeps=8,bond_dims=[32]*8,noises=[0.]*8,thrds=[1e-12]*8,tol=1e-12,iprint=0)
    else:ket=driver.load_mps('GS')
    identity=driver.get_identity_mpo()
    norm=float(driver.expectation(ket,identity,ket))
    sites=[0,1] if args.smoke else list(range(6))
    if not args.smoke:
        maps=json.loads((ROOT/'full_hamiltonian/index_maps.json').read_text())
        state_meta=json.loads((path/'energies.json').read_text())
        assert maps['cas20_local_to_full_zero_based'][:6]==maps['K6']
        sites=sorted(i for i,j in enumerate(state_meta['solver_site_to_local_zero_based']) if j<6)
    if args.smoke:
        from pyscf import ao2mo,fci
        g=ao2mo.restore(1,driver.g2e,driver.n_sites)
        _,ci=fci.direct_spin1.kernel(driver.h1e,g,driver.n_sites,driver.n_elec,ecore=driver.ecore,tol=1e-12)
        strs=fci.cistring.make_strings(range(driver.n_sites),driver.n_elec//2)
        j=int(np.argmax(abs(np.diag(np.asarray(ci)))**2))
        sites=[i for i in range(driver.n_sites) if int(strs[j])&(1<<i)]
        exact=float(abs(ci[j,j])**2)
        assert exact>.05,exact
    probability=float(driver.expectation(ket,projector_mpo(driver,sites),ket))/norm
    assert -1e-10<=probability<=1+1e-10,probability
    result=dict(sites_zero_based=sites,probability=probability,norm=norm,
        method='SU2 scalar MPO product of local n_alpha*n_beta projectors',
        interpretation='CAS20 probability of the specified fixed-double sector; not overlap with the reoptimized CAS14 state')
    result['fcidump_sha256']=hashlib.sha256((path/'FCIDUMP').read_bytes()).hexdigest()
    if not args.smoke:
        result['state_energy']=state_meta['energy_total']
        result['rdm1_sha256']=hashlib.sha256((path/'rdm1.npy').read_bytes()).hexdigest()
    if args.smoke:
        result.update(FCI_probability=exact,error=probability-exact)
        assert abs(probability-exact)<2e-5,result
        print('NONZERO_FCI_PROJECTOR_CHECK',probability,exact,flush=True)
        # Test the full 24-fermion six-projector expression at both endpoints.
        endpoints={}
        for ne,tag,occupations,expected in [(12,'SIXFULL',[2]*6,1.),(10,'SIXHOLE',[2]*5+[0],0.)]:
            driver.initialize_system(n_sites=6,n_elec=ne,spin=0,orb_sym=None)
            state=driver.get_random_mps(tag=tag,bond_dim=1,nroots=1,occs=occupations)
            print('SIX_PROJECTOR_ENDPOINT',tag,flush=True)
            ident=driver.get_identity_mpo();z=float(driver.expectation(state,ident,state))
            value=float(driver.expectation(state,projector_mpo(driver,list(range(6))),state))/z
            assert abs(value-expected)<1e-10,(tag,value)
            endpoints[tag]=value
        result['six_projector_endpoint_tests']=endpoints
    (path/'sector_probability.json').write_text(json.dumps(result,indent=2)+'\n')
    print('SECTOR_PROBABILITY',json.dumps(result))

if __name__=='__main__':main()
