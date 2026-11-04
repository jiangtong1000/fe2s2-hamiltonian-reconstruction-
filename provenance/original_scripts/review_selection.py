"""Audit stored chemical selection without selecting orbitals using energy."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from pyscf import scf

ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'orbitals'
    mol=scf.chkfile.load_mol(str(ROOT/'parent/scf.chk'))
    labels=mol.ao_labels(fmt=False)
    p=np.load(out/'character_meta_lowdin.npy')
    report=json.loads((out/'selection_report.json').read_text())
    assert (report['nfrozen'],report['ninactive'],report['nactive'],report['nvirtual'])==(44,34,20,564)
    def ids(atom,shell):
        return [i for i,x in enumerate(labels) if x[0]==atom and x[2]==shell]
    def spectrum(rows,cols):
        return np.linalg.eigvalsh(p[np.ix_(rows,cols)]@p[np.ix_(rows,cols)].T).tolist()
    core=[i for i,(_,el,shell,_) in enumerate(labels)
          if (el in ['Fe','S'] and shell in ['1s','2s','2p']) or (el=='C' and shell=='1s')]
    core_spectrum=spectrum(core,list(range(44)))
    semicore={str(a):{'frozen_weight':float(np.sum(p[ids(a,'3s')+ids(a,'3p'),:44]**2)),
                     'retained_projector_eigenvalues':spectrum(ids(a,'3s')+ids(a,'3p'),list(range(44,662)))} for a in [0,5]}
    bridge={str(a):spectrum(ids(a,'3p'),list(range(78,98))) for a in [1,4]}
    fe={str(a):spectrum(ids(a,'3d'),list(range(78,98))) for a in [0,5]}
    # Majority-character gates are declared diagnostics for this implementation,
    # not claimed historical thresholds. All directions must be represented.
    checks={'core_projector_min_above_0p99':min(core_spectrum)>.99,
            'semicore_retained_projector_min_above_0p99':all(min(x['retained_projector_eigenvalues'])>.99 for x in semicore.values()),
            'bridge_all_three_directions_above_0p5':all(min(x)>.5 for x in bridge.values()),
            'fe_all_five_directions_above_0p5':all(min(x)>.5 for x in fe.values()),
            'terminal_Fe_directed_weight_above_0p5':all(x['score']>.5 for x in report['ligand_assignments'] if x['terminal_fe_directed']),
            'complete_basis_orthonormal':report['orthonormality_max_error']<1e-10}
    rows=list(csv.DictReader((out/'orbital_table_initial.csv').open()))
    shell_groups=[(a,shell) for a in range(mol.natm) for shell in dict.fromkeys(x[2] for x in labels if x[0]==a)]
    for a,shell in shell_groups:
        weights=np.sum(p[ids(a,shell)]**2,axis=0)
        name=f'atom{a}_{mol.atom_symbol(a)}_{shell}_meta_lowdin_weight'
        for j,row in enumerate(rows): row[name]=float(weights[j])
    with (out/'orbital_character_shells.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    audit=dict(checks=checks,core_projector_eigenvalues=core_spectrum,Fe_semicore=semicore,
               bridge_3p_projector_eigenvalues=bridge,Fe_3d_projector_eigenvalues=fe,
               diagnostic_threshold_origin='Explicit implementation choices; not historical PM settings',
               selection_uses_energy=False,all_electrons=186,frozen_electrons=88,retained_electrons=98,
               retained_singlet_electrons=[49,49],CAS20_electrons=30,
               coefficients_sha256=hashlib.sha256((out/'mo_coeff_partition_initial.npy').read_bytes()).hexdigest())
    (out/'chemical_review.json').write_text(json.dumps(audit,indent=2)+'\n')
    assert all(checks.values()),audit
    report['status']='chemical_review_passed'
    report['chemical_review_artifact']='orbitals/chemical_review.json'
    (out/'selection_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__': main()
