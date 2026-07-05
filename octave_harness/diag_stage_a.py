"""Stage A/C alignment diagnostics (REPORT.md Priority 1).

Part 1 (Task 3): is process_cast.m step 8 (APPLY PITCH/ROLL CORRECTIONS)
a no-op on this cast? process_cast.m:317 gates it on length(p.tiltcor)>1;
default.m:309 defaults tiltcor to scalar 0 and the recorded LDEO p-struct
agrees (recorded_p_struct_attrs.txt:232) -- so we expect step07 == step08
bit-for-bit.

Part 2 (Task 4): row/column alignment diagnostics for the 14 m/s Stage A
residual -- see the functions below check_step8_noop().

Run from the LADCP repo root:  uv run python octave_harness/diag_stage_a.py
"""
from pathlib import Path
import sys

import numpy as np
import scipy.io as sio

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "octave_harness"))

DUMPS = REPO / "octave_harness" / "work" / "dumps"
DATA_DIR = REPO / "test_data" / "2015_P16N"


def _load(step: int):
    return sio.loadmat(
        DUMPS / f"step{step:02d}.mat", struct_as_record=False, squeeze_me=True
    )


def _f(x):
    return np.asarray(x, dtype=float)


def check_step8_noop() -> bool:
    s7, s8 = _load(7), _load(8)
    d7, d8, p8 = s7["d"], s8["d"], s8["p"]
    print(f"p.tiltcor = {p8.tiltcor!r}  (scalar => step 8 gate is False)")
    all_identical = True
    for name in ("ru", "rv", "rw"):
        a, b = _f(getattr(d7, name)), _f(getattr(d8, name))
        same = (np.isnan(a) & np.isnan(b)) | (a == b)
        ident = bool(same.all())
        all_identical &= ident
        print(f"d.{name}: shapes {a.shape}=={b.shape}  identical={ident}  "
              f"n_cells_differing={(~same).sum()}")
    print(f"\nSTEP-8 NO-OP: {all_identical}")
    return all_identical


if __name__ == "__main__":
    check_step8_noop()
