"""Small-system independent FCI checks of nesting and embedded expectations."""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))

import numpy as np
from pyscf import fci, lib
from integrals import (
    contract_fixed,
    active_eri,
    embedded_singlet_energy,
    determinant_energy,
)

lib.num_threads(1)
rng = np.random.default_rng(1287)
n = 6
h = rng.normal(size=(n, n))
h = (h + h.T) / 2
l = rng.normal(size=(9, n, n))
l = 0.08 * (l + l.transpose(0, 2, 1))
e0 = 1.23
f = [0]
i = [1]
a = [2, 3, 4, 5]
ret = list(range(1, n))
hr, er = contract_fixed(h, l, e0, f, ret)
lr = l[:, ret][:, :, ret]
ha, ea = contract_fixed(hr, lr, er, [0], [1, 2, 3, 4])
h_union, e_union = contract_fixed(h, l, e0, f + i, a)
assert abs(ha - h_union).max() < 1e-12
assert abs(ea - e_union) < 1e-12
eri = active_eri(l, a)
solver = fci.direct_spin0.FCI()
e, ci = solver.kernel(ha, eri, len(a), 2, ecore=ea)
dm1, dm2_pyscf = solver.make_rdm12(ci, len(a), 2)
dm2_block2 = dm2_pyscf.transpose(0, 2, 3, 1)
eembed = embedded_singlet_energy(h, l, e0, f + i, a, dm1, dm2_block2)
assert abs(e - eembed) < 1e-10, (e, eembed)
# Independent explicit full-ERI contraction for complex occupied coefficients.
c = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]
ca, cb = c[:, :3], c[:, 1:3]
da, db = ca @ ca.conj().T, cb @ cb.conj().T
d = da + db
g = np.einsum("xpq,xrs->pqrs", l, l)
ed = e0 + np.einsum("pq,qp", h, d) + 0.5 * np.einsum("pqrs,qp,sr", g, d, d)
ed -= 0.5 * sum(np.einsum("pqrs,sp,qr", g, ds, ds) for ds in [da, db])
ef = determinant_energy(h, l, e0, ca, cb)
assert abs(ed - ef) < 1e-10, (ed, ef)
print(
    "PASS sequential freezing; singlet FCI embedded cumulant energy; complex determinant energy"
)
