"""Transform full AO factors, contract only F, and derive both final CAS files.

Every retained cross-block coupling survives. Trial partition metadata is copied
beside the full Hamiltonian and is not interpreted as propagation constraints.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import numpy as np
from pyscf import lib,tools
from integrals import contract_fixed,active_eri

ROOT=Path(__file__).resolve().parents[1]
def dump_cas(path,h,l,e0,fixed,active,nelec):
    path.mkdir(exist_ok=True)
    ha,ea=contract_fixed(h,l,e0,fixed,active)
    eri=active_eri(l,active)
    np.save(path/'h1.npy',ha);np.save(path/'eri.npy',eri)
    (path/'ecore.txt').write_text(f'{ea:.16e}\n')
    tools.fcidump.from_integrals(str(path/'FCIDUMP'),ha,eri,len(active),nelec,
                                nuc=ea,ms=0,tol=1e-14,float_format=' %.16e')
    return dict(norb=len(active),nelec=nelec,ecore=ea,
                local_to_full_zero_based=list(active),fcidump_sha256=hashlib.sha256((path/'FCIDUMP').read_bytes()).hexdigest())
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    out=ROOT/'full_hamiltonian';out.mkdir(exist_ok=True)
    meta=json.loads((ROOT/'ao_factorization/metadata.json').read_text())
    maps=json.loads((ROOT/'orbitals/index_maps.json').read_text())
    c=np.load(ROOT/'orbitals/mo_coeff_full_final.npy')
    cf=np.load(ROOT/'orbitals/mo_coeff_frozen_core.npy')
    hcore=np.load(ROOT/'parent/hcore_ao.npy')
    parent=json.loads((ROOT/'parent/ao_metadata.json').read_text())
    s=np.load(ROOT/'parent/overlap_ao.npy')
    assert abs(c.T@s@c-np.eye(c.shape[1])).max()<1e-10
    assert abs(cf.T@s@c).max()<1e-10
    assert cf.shape[1]==44 and c.shape[1]==618
    la=np.load(ROOT/'ao_factorization/chol_ao_packed.npy',mmap_mode='r')
    rank,n=la.shape[0],c.shape[1]
    loose_rank=meta['ranks_at_threshold'].get('1e-08')
    if loose_rank is None: raise RuntimeError('Missing looser factorization checkpoint rank')
    l=np.lib.format.open_memmap(out/'chol.npy',mode='w+',dtype=np.float64,shape=(rank,n,n))
    h=c.T@hcore@c
    e0=float(parent['nuclear_repulsion']+2*np.trace(cf.T@hcore@cf))
    rows,cols=np.tril_indices(c.shape[0])
    start=0
    while start<rank:
        end=min(start+16,rank)
        if start<loose_rank<end:end=loose_rank
        a=np.empty((end-start,c.shape[0],c.shape[0]))
        a[:,rows,cols]=la[start:end];a[:,cols,rows]=la[start:end]
        ac=a@c;af=a@cf
        rr=c.T@ac;rf=c.T@af;ff=cf.T@af
        tr=np.trace(ff,axis1=1,axis2=2)
        h+=2*np.einsum('gpq,g->pq',rr,tr)-np.einsum('gpi,gqi->pq',rf,rf)
        e0+=float(2*np.dot(tr,tr)-np.einsum('gij,gij',ff,ff))
        assert np.max(abs(rr-rr.transpose(0,2,1)))<1e-10
        l[start:end]=rr
        start=end
        if end==loose_rank:
            np.save(out/'h1_loose.npy',h)
            (out/'h0_loose.txt').write_text(f'{e0:.16e}\n')
        if end%128==0:print('TRANSFORMED',end,'of',rank,flush=True)
    l.flush()
    np.save(out/'h1.npy',h)
    (out/'h0.txt').write_text(f'{e0:.16e}\n')
    (out/'nelec.txt').write_text('49 49\n')
    (out/'orbital_order.txt').write_text('[I | K6 | A14 | V]; explicit zero-based maps in index_maps.json\n')
    shutil.copy2(ROOT/'orbitals/index_maps.json',out/'index_maps.json')
    schema=dict(status='built_pending_validation',h1='physical, including F Coulomb/exchange; no XC or AFQMC shift',
        scalar='nuclear repulsion plus F exactly once',eri_convention='chemists',dtype='float64',
        energy_units='hartree',chol_layout=['nchol','norb','norb'],chol_shape=list(l.shape),
        ao_factorization_threshold=meta['threshold'],looser_threshold=1e-8,looser_rank=loose_rank,
        retained_electrons=98,retained_nalpha=49,retained_nbeta=49,
        frozen_core_orbitals=44,full_cross_block_couplings_retained=True,
        trial_fixed_occupations_restrict_propagation=False,
        backend='neutral; no AFQMC adapter',maps='index_maps.json')
    schema['cas20']=dump_cas(ROOT/'cas20',h,l,e0,maps['I'],maps['cas20_local_to_full_zero_based'],30)
    schema['cas14']=dump_cas(ROOT/'cas14',h,l,e0,maps['I']+maps['K6'],maps['A14'],18)
    # Same final NO basis and entire tight F contraction represented at looser rank.
    hl=np.load(out/'h1_loose.npy');el=float((out/'h0_loose.txt').read_text())
    schema['cas20_loose']=dump_cas(ROOT/'cas20_loose',hl,l[:loose_rank],el,maps['I'],maps['cas20_local_to_full_zero_based'],30)
    schema['cas14_loose']=dump_cas(ROOT/'cas14_loose',hl,l[:loose_rank],el,maps['I']+maps['K6'],maps['A14'],18)
    (out/'schema.json').write_text(json.dumps(schema,indent=2)+'\n')
    print('PRODUCTION_INTEGRALS_BUILT',json.dumps(schema))
if __name__=='__main__':main()
