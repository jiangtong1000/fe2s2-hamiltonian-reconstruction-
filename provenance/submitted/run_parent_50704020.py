"""High-spin physical parent and complete UKS-density natural orbital basis."""
from pathlib import Path
import hashlib
import json
import os
import platform
import numpy as np
from scipy.linalg import eigh
from pyscf import gto, dft, lib, __version__ as pyscf_version

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'parent'
def save(name, a): np.save(OUT/(name+'.npy'),a)
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    mol = gto.M(atom=str(ROOT/'geometry.xyz'),unit='Angstrom',basis='tzp-dkh',charge=-2,
                spin=10,symmetry=False,cart=False,verbose=4,max_memory=48000)
    assert mol.nelec == (98,88)
    mf = dft.UKS(mol).x2c()
    mf.xc='BP86'
    mf.grids.level=4
    mf.conv_tol=1e-10
    mf.conv_tol_grad=1e-6
    mf.max_cycle=160
    mf.chkfile=str(OUT/'scf.chk')
    energy=mf.kernel()
    if not mf.converged:
        raise RuntimeError('UKS not converged; checkpoint preserved for diagnosis')
    s=mf.get_ovlp()
    h=mf.get_hcore()
    dm=mf.make_rdm1()
    se,sv=eigh(s)
    if se.min()<1e-8:
        raise RuntimeError('AO overlap needs explicit linear-dependence policy before UNO construction')
    sh=(sv*np.sqrt(se))@sv.T
    si=(sv/np.sqrt(se))@sv.T
    no,u=eigh(sh@(dm[0]+dm[1])@sh)
    no,u=no[::-1],u[:,::-1]
    c=si@u
    phases=np.where(c[np.argmax(abs(c),axis=0),np.arange(c.shape[1])]<0,-1.,1.)
    c*=phases
    err=float(abs(c.T@s@c-np.eye(c.shape[1])).max())
    assert err<1e-10,err
    assert abs(no.sum()-186)<1e-8
    assert no.min()>-1e-8 and no.max()<2+1e-8
    for name,a in dict(mo_coeff=c,overlap_ao=s,hcore_ao=h,uks_mo_coeff=mf.mo_coeff,
                       uks_mo_occ=mf.mo_occ,uks_mo_energy=mf.mo_energy,
                       occupied_alpha=mf.mo_coeff[0][:,mf.mo_occ[0]>0],
                       occupied_beta=mf.mo_coeff[1][:,mf.mo_occ[1]>0],
                       density_alpha_ao=dm[0],density_beta_ao=dm[1],
                       parent_natural_occupations=no,
                       uks_alpha_to_uno=c.T@s@mf.mo_coeff[0],
                       uks_beta_to_uno=c.T@s@mf.mo_coeff[1]).items(): save(name,np.asarray(a))
    # Independent physical Hamiltonian expectation on the actual UKS determinant.
    j,k=mf.get_jk(mol,dm)
    e1=np.einsum('pq,qp',h,dm.sum(axis=0))
    e2=.5*np.einsum('pq,qp',j.sum(axis=0),dm.sum(axis=0))-.5*np.einsum('spq,sqp',k,dm)
    meta=dict(pyscf_version=pyscf_version,python=platform.python_version(),nao=mol.nao_nr(),
        nelec=list(mol.nelec),df_total_energy=float(energy),physical_determinant_total_energy=float(e1+e2+mol.energy_nuc()),
        nuclear_repulsion=float(mol.energy_nuc()),spin_square=mf.spin_square(),scf_converged=bool(mf.converged),
        xc=mf.xc,grid_level=mf.grids.level,conv_tol=mf.conv_tol,conv_tol_grad=mf.conv_tol_grad,
        overlap_smallest_eigenvalue=float(se.min()),uno_orthonormality_max_error=err,
        uno_source='spin-summed high-spin UKS 1-RDM',ao_functions_removed=0,
        sf_x2c=dict(class_name=type(mf.with_x2c).__name__,xuncontract=mf.with_x2c.xuncontract,
                    approximation=mf.with_x2c.approx,two_electron_picture_change=False),
        basis=mol._basis,ao_labels=mol.ao_labels(),spherical=True,
        geometry_sha256=hashlib.sha256((ROOT/'geometry.xyz').read_bytes()).hexdigest())
    (OUT/'ao_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    np.savetxt(OUT/'parent_natural_occupations.txt',np.column_stack([np.arange(len(no)),no]),fmt=['%d','%.14f'])
    print('PARENT_COMPLETE',json.dumps({k:meta[k] for k in ['nao','df_total_energy','physical_determinant_total_energy','uno_orthonormality_max_error']}))
if __name__=='__main__': main()
