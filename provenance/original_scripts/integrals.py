"""Backend-neutral chemists' ERI contractions; symmetric unshifted L[g,p,q]."""
import numpy as np

def batches(l, size=32):
    for start in range(0,len(l),size):
        yield np.asarray(l[start:start+size])

def contract_fixed(h,l,e0,fixed,active):
    """Exact fixed-double/vacuum subspace contraction of given factors."""
    fixed=np.asarray(fixed,dtype=int)
    active=np.asarray(active,dtype=int)
    assert len(set(fixed)&set(active))==0
    ha=h[np.ix_(active,active)].copy()
    ea=float(e0+2*np.trace(h[np.ix_(fixed,fixed)]))
    for b in batches(l):
        kk=b[:,fixed][:,:,fixed]
        ak=b[:,active][:,:,fixed]
        aa=b[:,active][:,:,active]
        tr=np.trace(kk,axis1=1,axis2=2)
        ha+=2*np.einsum('gpq,g->pq',aa,tr)-np.einsum('gpi,gqi->pq',ak,ak)
        ea+=float(2*np.dot(tr,tr)-np.einsum('gij,gij',kk,kk))
    return ha,ea

def active_eri(l,active):
    a=np.asarray(active,dtype=int)
    eri=np.zeros((len(a),)*4)
    for b in batches(l):
        aa=b[:,a][:,:,a]
        eri+=np.einsum('gpq,grs->pqrs',aa,aa,optimize=True)
    return eri

def determinant_energy(h,l,e0,ca,cb):
    """Physical expectation for real or complex occupied spin coefficients."""
    da=ca@ca.conj().T
    db=cb@cb.conj().T
    d=da+db
    energy=complex(e0+np.einsum('pq,qp',h,d))
    for b in batches(l):
        occupied_blocks=[cs.conj().T@(b@cs) for cs in (ca,cb)]
        trace=sum(np.trace(x,axis1=1,axis2=2) for x in occupied_blocks)
        energy+=.5*np.dot(trace,trace)
        for block in occupied_blocks:
            energy-=.5*np.einsum('gij,gji',block,block)
    assert abs(energy.imag)<1e-7,energy
    return float(energy.real)

def embedded_singlet_energy(h,l,e0,fixed,active,rdm1,rdm2):
    """Full-density plus active-cumulant route independent of contract_fixed.

    dm2[p,r,s,q] = sum_st <p_s^+ r_t^+ s_t q_s> (block2 convention).
    """
    a=np.asarray(active,dtype=int)
    k=np.asarray(fixed,dtype=int)
    occupied=np.concatenate([k,a])
    d=np.zeros((len(occupied),len(occupied)))
    d[np.arange(len(k)),np.arange(len(k))]=2
    d[len(k):,len(k):]=rdm1
    e=float(e0+np.einsum('pq,qp',h[np.ix_(occupied,occupied)],d))
    cumulant=rdm2.transpose(0,3,1,2)-np.einsum('pq,rs->pqrs',rdm1,rdm1)+.5*np.einsum('ps,rq->pqrs',rdm1,rdm1)
    for b in batches(l):
        bo=b[:,occupied][:,:,occupied]
        trace=np.einsum('gpq,qp->g',bo,d)
        ld=bo@d
        e+=.5*np.dot(trace,trace)-.25*np.einsum('gpq,gqp',ld,ld)
        aa=b[:,a][:,:,a]
        e+=.5*np.einsum('gpq,grs,pqrs',aa,aa,cumulant,optimize=True)
    return float(e)

def direct_cas(mol,hcore,c,core,active):
    """Independent direct AO integrals and AO J/K contraction for small CAS."""
    from pyscf import scf, ao2mo
    core=np.asarray(core,dtype=int)
    active=np.asarray(active,dtype=int)
    cc=c[:,core]
    ca=c[:,active]
    d=2*cc@cc.T
    j,k=scf.hf.get_jk(mol,d)
    v=j-.5*k
    h=ca.T@(hcore+v)@ca
    e0=float(mol.energy_nuc()+np.einsum('pq,qp',hcore+.5*v,d))
    eri=ao2mo.restore(1,ao2mo.kernel(mol,ca),len(active))
    return h,eri,e0
