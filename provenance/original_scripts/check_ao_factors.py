"""Bounded-memory direct molecular AO quartet checks on serialized factors."""
from pathlib import Path
import json
import os
import numpy as np
from pyscf import scf,lib

ROOT=Path(__file__).resolve().parents[1]
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    meta=json.loads((ROOT/'ao_factorization/metadata.json').read_text())
    l=np.load(ROOT/'ao_factorization/chol_ao_packed.npy',mmap_mode='r')
    mol=scf.chkfile.load_mol(str(ROOT/'parent/scf.chk'))
    offsets=mol.ao_loc_nr();n=mol.nao_nr();rng=np.random.default_rng(20261005)
    quartets=[]
    for _ in range(32):
        p,q,r,s=rng.integers(0,n,4).tolist()
        quartets += [(p,q,p,q),(p,p,q,q),(p,q,r,s)]
    def pair(p,q):
        a,b=max(p,q),min(p,q)
        return a*(a+1)//2+b
    pair_ids=sorted({pair(p,q) for quartet in quartets for p,q in [quartet[:2],quartet[2:]]})
    # Stream row blocks once; fancy indexing a network-backed memmap by columns
    # creates many small, strided reads. Keep only the small selected panel.
    panel=np.empty((len(l),len(pair_ids)))
    shape,dtype,offset=l.shape,l.dtype,l.offset
    l._mmap.close()
    with (ROOT/'ao_factorization/chol_ao_packed.npy').open('rb') as f:
        f.seek(offset)
        for start in range(0,shape[0],32):
            count=min(32,shape[0]-start)
            raw=f.read(count*shape[1]*dtype.itemsize)
            block=np.frombuffer(raw,dtype=dtype).reshape(count,shape[1])
            panel[start:start+count]=block[:,pair_ids]
    lookup={p:i for i,p in enumerate(pair_ids)}
    samples=[]
    for k,(p,q,r,s) in enumerate(quartets):
        indices=[p,q,r,s]
        shells=[int(np.searchsorted(offsets,i,side='right')-1) for i in indices]
        local=tuple(i-offsets[sh] for i,sh in zip(indices,shells))
        exact=float(mol.intor_by_shell('int2e_sph',tuple(shells))[local])
        recovered=float(panel[:,lookup[pair(p,q)]]@panel[:,lookup[pair(r,s)]])
        samples.append(dict(indices=indices,family=['exchange_pair_diagonal','Coulomb','general'][k%3],
            direct=exact,factorized=recovered,error=recovered-exact))
    errors=np.array([x['error'] for x in samples]);tol=meta['threshold']*1.05+1e-11
    result=dict(seed=20261005,sampling_rule='32 uniform AO quartets, each augmented with exchange-pair diagonal and Coulomb patterns',
        maximum_error=float(abs(errors).max()),rms_error=float(np.sqrt(np.mean(errors**2))),
        tolerance=tol,passed=bool(abs(errors).max()<tol),samples=samples,
        caveat='AO sample check only; full MO virtual couplings and energy checks remain separate requirements')
    (ROOT/'validation/ao_factor_samples.json').write_text(json.dumps(result,indent=2)+'\n')
    assert result['passed'],result
    print('AO_FACTOR_SAMPLES_PASSED',result['maximum_error'])

if __name__=='__main__':main()
