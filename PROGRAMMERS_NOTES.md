# Programmer's Notes

Technical reference for developers working on this codebase. Read alongside the MATLAB source in `docs/legacy/` — the Python implementation is designed to be traceable to the MATLAB reference line-by-line.

## Architecture

Five layers in dependency order. Only the first is implemented.

```
Ingestion  ──▶  Transforms  ──▶  Solution  ──▶  QA / Diagnostics  ──▶  CLI / API
(done)          (stub)           (stub)          (stub)                 (stubs)
```

### Why two-layer ingestion

`src/ladcp/ingestion/` contains two internal modules:

- `_pd0.py` — low-level binary parser. Input: `bytes`. Output: `list[dict]`, one dict per ensemble, with raw-typed values (int16 counts, not physical units) where that matches the binary format, and physical units only where the conversion is part of the format spec (e.g. heading is already in 0.01° units from the instrument).
- `rdi.py` — public API. Calls `_pd0.parse_pd0()`, assembles dict-lists into numpy arrays, applies the 0.001 m/s velocity scale and NaN substitution, and returns `RDIData`.

This separation keeps the byte-parsing logic independently testable with synthetic byte buffers, while `rdi.py` tests can use minimal parsed-dict fixtures. Integration tests run against real `.000` files gated by `TEST_DATA_DIR`.

## Key Design Decisions

### Julian day convention: midnight-based (matching `julian.m`)

The MATLAB reference uses a non-standard Julian day convention where JD starts at **midnight**, not astronomical noon. This is implemented in `_pd0._to_julian()` using the Fliegel/Van Flandern algorithm:

```python
j = (146097 * c) // 4 + (1461 * yr) // 4 + (153 * mo + 2) // 5 + day + 1721119
return float(j) + hour_frac / 24.0
```

The equivalent MATLAB `docs/legacy/julian.m` key line is `j = j + h/24`, explicitly not `j + (h-12)/24` (which would be the Meeus astronomical noon-based formula). The difference is exactly 0.5 JD — using the wrong algorithm shifts every timestamp by 12 hours, corrupting DL/UL clock-drift corrections downstream.

The unit test assertion uses tight tolerance: `abs(result - 2458428.5) < 1e-4`.

### Byte offsets in the PD0 format

The Teledyne RDI PD0 format is documented in `docs/legacy/loadrdi.m`. Field positions inside each block are *relative offsets* — every block starts at a base address found in the offset table at bytes 6+ of the ensemble header. The offsets below are relative to the block start (byte 0 = type ID byte 0):

**Fixed leader (type 0x0000), `rdflead()` in loadrdi.m:**
- `nbin` at +7 (1 byte)
- `npng` at +8, `blen_cm` at +10, `blnk_cm` at +12
- `dist_cm` (distance to first bin centre) at +30 (2 bytes) — NOT +32
- `serial` at +40 — NOT +42

**Variable leader (type 0x0080), `rdvlead()` in loadrdi.m:**
- Timestamp (7 bytes: year, month, day, hour, min, sec, hundredths) at +2
- After timestamp: 3-byte skip (not 5-byte), then `sound_vel` at +12 — NOT +14
- `heading` at +16, pitch at +18, roll at +20 — NOT +18/20/22
- `salinity` at +22, `temp` at +24

These offsets were verified against the MATLAB reference by tracing `rdflead()`/`rdvlead()` in loadrdi.m. An earlier version of the implementation plan had transcription errors at +32/+42 and +14/+18 — if you see those values in any planning document, the corrected values above govern.

### Ensemble resync on bad length

When `parse_pd0()` encounters an ensemble with a declared byte-count that would read past the end of the buffer, it advances one byte and retries (`offset += 1; continue`). This handles:
- Partial ensembles at file boundaries
- Bit errors in the length word
- Files concatenated with padding bytes

An earlier `break` on bad length would silently truncate all remaining ensembles in the file.

### Coordinate frame assumption

`RDIData.u/v/w/e` are labelled as Earth-frame (East/North/Up/Error) only when the instrument was configured to output Earth-frame velocity — Teledyne RDI EX command `EX=11xxx`. If the instrument was configured for beam-frame or instrument-frame output, the arrays still load but `u` and `v` contain beam or instrument coordinates. The parser does not check the EX byte. This is documented in `_types.py`.

### Velocity scaling and NaN

Velocity fields in PD0 are int16, unit 0.001 m/s, with the sentinel value -32768 meaning "bad data." The parser converts -32768 to `np.nan` and scales everything else by 0.001. Bottom-track ranges are uint16 at 0.01 m/LSB with sentinel 0.

## Module Map

```
src/ladcp/
├── __init__.py              version = "0.1.0"
├── cli.py                   Click app: `ladcp process` and `ladcp check` (stubs)
├── ingestion/
│   ├── __init__.py          exports load_rdi
│   ├── _pd0.py              parse_pd0(), internal helpers
│   ├── _types.py            RDIData dataclass
│   └── rdi.py               load_rdi(path) → RDIData
├── transforms/
│   └── beam2earth.py        janus5beam2earth() stub (NotImplementedError)
├── solution/
│   └── shear.py             compute_shear() → ShearProfile
└── qa/
    └── diagnostics.py       tilt_heading_plot() stub
```

## Testing Approach

Tests live in `tests/`. The structure mirrors three levels of confidence:

**Unit tests** (`tests/test_pd0_parser.py`) — synthetic byte buffers built with helper functions (`_make_minimal_ensemble`, etc.). No external files needed. Run in CI unconditionally.

**Integration tests** (`tests/integration/`) — real PD0 files. Skipped when files are absent. Gated by `TEST_DATA_DIR` env var pointing to the `test_data/` directory. Two suites:
- `test_pd0_cast002.py` — I7N 2018 cast 002. Raw files not yet available; all 8 tests skip.
- `test_pd0_p16n_cast003.py` — 2015 P16N cast 003. Raw files available at `test_data/2015_P16N/`. 9 tests pass.

Run with: `TEST_DATA_DIR=test_data uv run pytest`

The `conftest.py` `test_data_dir` fixture provides the base path. Each integration test also checks for its specific file and calls `pytest.skip()` if absent — doubly safe.

## Extending the Ingestion Layer

To add support for a new PD0 block type:

1. Add a new type-ID dispatch case in `_pd0.parse_pd0()` (look for the `if type_id ==` chain).
2. Write a `_read_<name>()` helper following the same `data[offset:offset+N]` pattern.
3. Add the field to `RDIData` in `_types.py` with shape documentation.
4. Assemble the array in `rdi.load_rdi()`.
5. Add a unit test with a synthetic buffer and an integration assertion.

## Implementing the Transform Layer

The authoritative reference is `docs/legacy/ADCPtools/janus5beam2earth.m`. The Python stub signature in `transforms/beam2earth.py` should match the MATLAB function's parameters. Key options to replicate: `Gimbaled` (use pitch/roll from variable leader) and `Binmap` (correct for beam angle and tilt).

## Running Linter

```bash
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
```

Rules: E, F, I (imports), NPY (numpy), UP (pyupgrade). Line length is ruff's default (88).
