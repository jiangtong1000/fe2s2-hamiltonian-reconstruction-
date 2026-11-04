"""Independent molecular checks, including external-virtual couplings.

Run inside an allocation after production Hamiltonian creation. CAS state
expectations are included when final RDMs exist, and reported pending otherwise.
"""
from pathlib import Path
import sys
import json
import os
import numpy as np
from pyscf import scf,ao2mo,tools,lib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'build'))
from integrals import direct_cas,contract_fixed,active_eri,determinant_energy,embedded_singlet_energy
from cross_basis_energy import evaluate as cross_basis_energy

ROOT=Path(__file__).resolve().parents[1]
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    mol=scf.chkfile.load_mol(str(ROOT/'parent/scf.chk'))
    mol.max_memory=float(os.environ.get('VALIDATION_MEMORY_MB','24000'))
    s=np.load(ROOT/'parent/overlap_ao.npy')
    hao=np.load(ROOT/'parent/hcore_ao.npy')
    c=np.load(ROOT/'orbitals/mo_coeff_full_final.npy')
    cf=np.load(ROOT/'orbitals/mo_coeff_frozen_core.npy')
    c_all=np.column_stack([cf,c])
    maps=json.loads((ROOT/'full_hamiltonian/index_maps.json').read_text())
    h=np.load(ROOT/'full_hamiltonian/h1.npy')
    l=np.load(ROOT/'full_hamiltonian/chol.npy',mmap_mode='r')
    e0=float((ROOT/'full_hamiltonian/h0.txt').read_text())
    schema=json.loads((ROOT/'full_hamiltonian/schema.json').read_text())
    rng=np.random.default_rng(20261005)
    metrics={'status':'in_progress','tolerances':{'orthonormality':1e-10,'internal_integral':1e-10,
        'direct_energy_Eh':1e-5,'sampled_direct_ERI_Eh':1e-6,'embedded_rdm_Eh':1e-7},'checks':{}}
    checks=metrics['checks']
    def record(name,data):
        checks[name]=data
        (ROOT/'validation/validation_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    err=float(abs(c_all.T@s@c_all-np.eye(c_all.shape[1])).max())
    record('complete_orthonormality',dict(max_error=err,passed=err<1e-10))
    transformations=np.load(ROOT/'orbitals/transformations.npz')
    rebuilt=np.load(ROOT/'parent/mo_coeff.npy')@transformations['selection_from_parent']@transformations['final_partition_from_initial']
    transform_error=float(abs(rebuilt-c_all).max())
    record('saved_basis_transformations',dict(maximum_coefficient_error=transform_error,
        convention='C_final_complete = C_parent_UNO @ selection_from_parent @ final_partition_from_initial',
        passed=transform_error<1e-10))
    blocks=[maps[k] for k in ['I','K6','A14','V']]
    assert sorted(sum(blocks,[]))==list(range(c.shape[1]))
    assert maps['nelec_full']==2*len(maps['I'])+30==2*(len(maps['I'])+6)+18
    assert abs(c[:,maps['cas20_local_to_full_zero_based']]-np.load(ROOT/'orbitals/mo_coeff_20_final.npy')).max()<1e-12
    assert abs(c[:,maps['A14']]-np.load(ROOT/'orbitals/mo_coeff_14.npy')).max()<1e-12
    record('partition_and_nesting',dict(passed=True,nfrozen=cf.shape[1],nretained=c.shape[1]))
    cross_basis=cross_basis_energy()
    record('cas20_cross_basis_trial',cross_basis)
    if 'production_within_1e_5_of_saved_trial' in cross_basis:
        record('cas20_variational_convergence_crosscheck',dict(
            energy_error_lower_bound=cross_basis['production_energy_error_lower_bound'],
            passed=cross_basis['production_within_1e_5_of_saved_trial'],
            note='A necessary convergence check against an existing trial in this same Hamiltonian; passing alone does not certify the ground state'))
    cas_data={}
    for name,core,active,ne in [('cas20',maps['I'],maps['cas20_local_to_full_zero_based'],30),
                                ('cas14',maps['I']+maps['K6'],maps['A14'],18)]:
        path=ROOT/name
        parsed=tools.fcidump.read(str(path/'FCIDUMP'))
        assert parsed['NORB']==len(active) and parsed['NELEC']==ne and parsed['MS2']==0
        hh,ee=contract_fixed(h,l,e0,core,active)
        gg=active_eri(l,active)
        gg_text=ao2mo.restore(1,parsed['H2'],len(active))
        internal=dict(h1_max_error=float(abs(hh-parsed['H1']).max()),
                      eri_max_error=float(abs(gg-gg_text).max()),scalar_error=float(abs(ee-parsed['ECORE'])))
        internal['passed']=max(internal.values())<1e-10
        record(name+'_internal_integrals',internal)
        # Independent direct AO integral and inactive J/K route, including F once.
        direct_core=list(range(cf.shape[1]))+[cf.shape[1]+i for i in core]
        direct_active=[cf.shape[1]+i for i in active]
        hd,gd,ed=direct_cas(mol,hao,c_all,direct_core,direct_active)
        np.savez(ROOT/'validation'/(name+'_direct_integrals.npz'),h1=hd,eri=gd,ecore=ed)
        direct_errors=dict(h1_max_error=float(abs(hd-hh).max()),
                eri_max_error=float(abs(gd-gg).max()),scalar_error=float(ed-ee))
        direct_errors['passed']=max(direct_errors['h1_max_error'],direct_errors['eri_max_error'])<1e-6 and abs(ed-ee)<1e-5
        record(name+'_direct_integrals',direct_errors)
        cas_data[name]=(hh,gg,ee)
        if (path/'rdm2.npy').exists() and (path/'energies.json').exists():
            stages=json.loads((path/'convergence.json').read_text())
            energy_change=float(stages[-1]['energy_total']-stages[-2]['energy_total'])
            record(name+'_solver_convergence',dict(last_bond_stage_energy_change=energy_change,
                target_Eh=1e-5,passed=abs(energy_change)<1e-5,
                caveat='Successive bond-stage energy change is a convergence indicator, not a rigorous absolute-error bound'))
            dm1=np.load(path/'rdm1.npy');dm2=np.load(path/'rdm2.npy')
            energy=float(json.loads((path/'energies.json').read_text())['energy_total'])
            embedded=embedded_singlet_energy(h,l,e0,core,active,dm1,dm2)
            direct=float(ed+np.einsum('pq,pq',hd,dm1)+.5*np.einsum('pqrs,prsq',gd,dm2))
            record(name+'_embedded_and_direct_energies',dict(energy_solver=energy,energy_embedded=embedded,
                energy_direct=direct,embedded_error=embedded-energy,direct_error=direct-energy,
                passed=abs(embedded-energy)<1e-7 and abs(direct-energy)<1e-5))
            loose_path=ROOT/(name+'_loose')
            loose_h=np.load(loose_path/'h1.npy');loose_g=np.load(loose_path/'eri.npy')
            loose_e=float((loose_path/'ecore.txt').read_text())
            tight_expectation=float(ee+np.einsum('pq,pq',hh,dm1)+.5*np.einsum('pqrs,prsq',gg,dm2))
            loose_expectation=float(loose_e+np.einsum('pq,pq',loose_h,dm1)+.5*np.einsum('pqrs,prsq',loose_g,dm2))
            sensitivity=loose_expectation-tight_expectation
            record(name+'_factorization_sensitivity',dict(tight_state_tight_energy=tight_expectation,
                tight_state_loose_energy=loose_expectation,loose_minus_tight=sensitivity,
                note='Fixed-state expectation sensitivity; not a bound on relaxed-state or solver error',
                passed=abs(sensitivity)<1e-5))
        else:record(name+'_embedded_and_direct_energies',dict(status='pending final RDMs'))
    # Explicit small four-index algebra, independent of factor contraction helper.
    h20,g20,e20=cas_data['cas20'];keep=np.arange(6,20)
    h14=h20[np.ix_(keep,keep)].copy();e14=float(e20+2*np.trace(h20[:6,:6]))
    for i in range(6):
        h14+=2*g20[np.ix_(keep,keep,[i],[i])][:,:,0,0]-g20[np.ix_(keep,[i],[i],keep)][:,0,0,:]
        for j in range(6):e14+=2*g20[i,i,j,j]-g20[i,j,j,i]
    h14u,g14u,e14u=cas_data['cas14']
    seq=dict(h1_max_error=float(abs(h14-h14u).max()),scalar_error=float(abs(e14-e14u)),
             eri_max_error=float(abs(g20[6:,6:,6:,6:]-g14u).max()))
    seq['passed']=max(seq.values())<1e-10
    record('sequential_freezing',seq)
    # Full physical h1 after F contraction has an independent direct AO check.
    df=2*cf@cf.T;jf,kf=scf.hf.get_jk(mol,df);vf=jf-.5*kf
    h_direct=c.T@(hao+vf)@c
    e0_direct=float(mol.energy_nuc()+np.einsum('pq,qp',hao+.5*vf,df))
    one_body_error=float(abs(h-h_direct).max())
    record('full_physical_one_body',dict(h1_max_error=one_body_error,scalar_error=e0-e0_direct,
        passed=one_body_error<1e-6 and abs(e0-e0_direct)<1e-5))
    # Stratified external-virtual pool and direct AO-to-small-pool transformation.
    pool=[]; V=maps['V']; external_report={}
    for block in ['I','K6','A14']:
        b=maps[block]
        hv=h[np.ix_(b,V)]
        ib,iv=np.unravel_index(abs(hv).argmax(),hv.shape)
        picks=[b[ib],int(rng.choice(b))]
        vpicks=[V[iv],int(rng.choice(V))]
        pool+=picks+vpicks
        external_report[block]=dict(h1_frobenius_norm=float(np.linalg.norm(hv)),
            h1_max_abs=float(abs(hv).max()),h1_example_indices=[b[ib],V[iv]],
            h1_direct_error_max=float(abs(hv-h_direct[np.ix_(b,V)]).max()))
    pool=list(dict.fromkeys(pool))
    cp=c[:,pool]
    g_pool=ao2mo.restore(1,ao2mo.kernel(mol,cp),len(pool))
    reconstructed=active_eri(l,pool)
    samples=[]
    for block in ['I','K6','A14']:
        bs=[pool.index(i) for i in pool if i in maps[block]]
        vs=[pool.index(i) for i in pool if i in V]
        for family in ['paqb','pqab','paqr']:
            ids=[]
            for p in bs:
                for q in bs:
                    for a in vs:
                        for b in (bs if family=='paqr' else vs):
                            ids.append((p,a,q,b) if family in ['paqb','paqr'] else (p,q,a,b))
            values=np.array([g_pool[ix] for ix in ids])
            errs=np.array([reconstructed[ix]-g_pool[ix] for ix in ids])
            largest=int(abs(values).argmax());ix=ids[largest]
            samples.append(dict(block=block,family=family,count=len(ids),
                max_abs_direct=float(abs(values).max()),rms_direct=float(np.sqrt(np.mean(values**2))),
                max_error=float(abs(errs).max()),rms_error=float(np.sqrt(np.mean(errs**2))),
                representative_full_indices=[pool[i] for i in ix],direct_value=float(values[largest]),
                reconstructed_value=float(reconstructed[ix])))
    coupling_pass=all(x['h1_max_abs']>1e-10 and x['h1_direct_error_max']<1e-6 for x in external_report.values())
    coupling_pass &= all(x['max_abs_direct']>1e-10 and x['max_error']<1e-6 for x in samples)
    record('external_virtual_couplings',dict(one_body_blocks=external_report,eri_families=samples,
        pool_full_indices=pool,seed=20261005,sampling_rule='largest physical h1 entry and one random occupied/active column per block; associated plus random V columns; all listed families within this bounded pool',
        arrays_reloaded_from_disk=True,nonzero_threshold=1e-10,passed=bool(coupling_pass),
        caveat='sampled errors are not global bounds; failed zero blocks require physical diagnosis, not fabricated couplings'))
    # Transform the actual UKS occupied spaces, then explicitly construct an
    # F-compatible determinant. The parent spin blocks need not doubly occupy F.
    parent=json.loads((ROOT/'parent/ao_metadata.json').read_text())
    compatible=[];parent_transforms={};projection={}
    for spin,ne in [('alpha',98),('beta',88)]:
        occupied=np.load(ROOT/'parent'/f'occupied_{spin}.npy')
        transformed=c_all.T@s@occupied
        retained=transformed[cf.shape[1]:]
        uu,singular,_=np.linalg.svd(retained,full_matrices=False)
        q=uu[:,:ne-cf.shape[1]]
        compatible.append(q)
        parent_transforms[spin+'_actual_occupied_in_complete_final_basis']=transformed
        parent_transforms[spin+'_compatible_occupied_in_retained_basis']=q
        projection[spin]=dict(F_occupation=float(np.sum(transformed[:cf.shape[1]]**2)),
            F_holes=float(cf.shape[1]-np.sum(transformed[:cf.shape[1]]**2)),
            retained_singular_values=singular.tolist(),retained_occupied_count=q.shape[1])
    np.savez(ROOT/'validation/parent_determinant_transforms.npz',**parent_transforms)
    ca,cb=[np.column_stack([cf,c@q]) for q in compatible]
    assert ca.shape[1]==98 and cb.shape[1]==88
    assert abs(ca.T@s@ca-np.eye(98)).max()<1e-10
    assert abs(cb.T@s@cb-np.eye(88)).max()<1e-10
    dm=np.array([ca@ca.T,cb@cb.T]);j,k=scf.hf.get_jk(mol,dm)
    compatible_ao=float(mol.energy_nuc()+np.einsum('pq,qp',hao,dm.sum(axis=0))+.5*np.einsum('pq,qp',j.sum(axis=0),dm.sum(axis=0))-.5*np.einsum('spq,sqp',k,dm))
    compatible_mo=determinant_energy(h,l,e0,*compatible)
    record('parent_compatible_determinant',dict(projection=projection,
        construction='F doubly occupied; retain dominant Nspin-nF left singular vectors of projected parent occupied coefficients',
        compatible_AO_energy=compatible_ao,compatible_MO_energy=compatible_mo,
        unprojected_parent_physical_energy=parent['physical_determinant_total_energy'],
        compatible_minus_unprojected_parent=compatible_ao-parent['physical_determinant_total_energy'],
        AO_MO_error=compatible_mo-compatible_ao,passed=abs(compatible_mo-compatible_ao)<1e-5,
        note='High-spin integral validation only; neither the DFT total energy nor a singlet trial'))
    # Determinants with external mixing test information invisible to CAS recovery.
    hl=np.load(ROOT/'full_hamiltonian/h1_loose.npy')
    el=float((ROOT/'full_hamiltonian/h0_loose.txt').read_text())
    determinants=[]
    for test in range(3):
        occa=np.eye(c.shape[1])[:,:49].copy();occb=occa.copy()
        virtual_choices=rng.choice(V,size=3,replace=False)
        for b,a in zip(['I','K6','A14'],virtual_choices):
            candidates=[i for i in maps[b] if i<49]
            i=int(rng.choice(candidates));a=int(a)
            angle=.15*(test+1)
            occa[i,i]=np.cos(angle);occa[a,i]=np.sin(angle)
            occb[i,i]=np.cos(angle/2);occb[a,i]=-np.sin(angle/2)
        ca=np.column_stack([cf,c@occa]);cb=np.column_stack([cf,c@occb])
        assert abs(ca.T@s@ca-np.eye(93)).max()<1e-10
        assert abs(cb.T@s@cb-np.eye(93)).max()<1e-10
        da=ca@ca.T;db=cb@cb.T;dm=np.array([da,db]);j,k=scf.hf.get_jk(mol,dm)
        direct=float(mol.energy_nuc()+np.einsum('pq,qp',hao,da+db)+.5*np.einsum('pq,qp',j.sum(axis=0),da+db)-.5*np.einsum('spq,sqp',k,dm))
        tight=determinant_energy(h,l,e0,occa,occb)
        loose=determinant_energy(hl,l[:schema['looser_rank']],el,occa,occb)
        determinants.append(dict(test=test,ao_direct=direct,mo_tight=tight,mo_loose=loose,
                                 tight_error=tight-direct,loose_to_tight=loose-tight,
                                 passed=abs(tight-direct)<1e-5 and abs(loose-tight)<1e-5))
    record('virtual_mixed_determinants',dict(tests=determinants,passed=all(x['passed'] for x in determinants)))
    pending=any(v.get('status','').startswith('pending') for v in checks.values())
    fails=[k for k,v in checks.items() if v.get('passed') is False]
    metrics.update(status='failed' if fails else 'pending_states' if pending else 'numerical_checks_passed',failed_checks=fails)
    (ROOT/'validation/validation_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print('VALIDATION_STATUS',metrics['status'],fails)
if __name__=='__main__':main()
