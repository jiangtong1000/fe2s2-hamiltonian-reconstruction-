"""Direct-integral CAS20 for selection/DMRG; not the production source."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from pyscf import scf,lib,tools
from integrals import direct_cas

ROOT=Path(__file__).resolve().parents[1]
def main():
    lib.num_threads(int(os.environ.get('OMP_NUM_THREADS','1')))
    report=json.loads((ROOT/'orbitals/selection_report.json').read_text())
    if report['status']!='chemical_review_passed':
        raise RuntimeError('Inspect selection diagnostics before accepting the chemical CAS')
    out=ROOT/'cas20_provisional'; out.mkdir(exist_ok=True)
    mol=scf.chkfile.load_mol(str(ROOT/'parent/scf.chk'))
    mol.max_memory=int(os.environ.get('CAS_BUILD_MEMORY_MB','6000'))
    c=np.load(ROOT/'orbitals/mo_coeff_partition_initial.npy')
    hcore=np.load(ROOT/'parent/hcore_ao.npy')
    h,eri,e0=direct_cas(mol,hcore,c,list(range(78)),list(range(78,98)))
    np.save(out/'h1.npy',h); np.save(out/'eri.npy',eri)
    tools.fcidump.from_integrals(str(out/'FCIDUMP'),h,eri,20,30,nuc=e0,ms=0,
                                tol=1e-14,float_format=' %.16e')
    occ=np.full(20,2); occ[[2,3,4,5,6,13,14,15,16,17]]=1
    np.save(out/'initial_occupations.npy',occ)
    meta=dict(stage='provisional_direct_AO',norb=20,nelec=30,ecore=float(e0),
              selection_source='orbitals/selection_report.json',
              fcidump_sha256=hashlib.sha256((out/'FCIDUMP').read_bytes()).hexdigest(),
              basis_artifact='orbitals/mo_coeff_20_initial.npy')
    (out/'hamiltonian_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('PROVISIONAL_CAS_READY',json.dumps(meta))
if __name__=='__main__': main()
