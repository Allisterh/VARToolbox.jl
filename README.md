# VARToolbox.jl

A Julia port of **Ambrogio Cesa-Bianchi's VAR Toolbox** (MATLAB, v3.1):
reduced-form VAR estimation, five identification schemes, impulse responses,
forecast error variance decompositions, historical decompositions, bootstrap
error bands, local projections and ARDL models.

Licensed **GPL-3.0**, inherited from the original. See
[`PROVENANCE.md`](PROVENANCE.md) for the upstream commit this derives from.

## Status

The estimation and inference pipeline is complete and validated against the
MATLAB original. **1,602 tests pass**, almost all of them exact
golden-master comparisons at `rtol = 1e-10` against output dumped from
MATLAB R2026a.

All four published replications shipped with the toolbox reproduce exactly:

| replication | identification | agreement |
|---|---|---|
| Stock & Watson (2001, JEP) | recursive | `1e-10` |
| Blanchard & Quah (1989, AER) | long-run | `1e-10` |
| Uhlig (2005, JME) | sign restrictions | `1e-10` |
| Gertler & Karadi (2015, AEJ) | external instrument | `1e-6` [^cond] |

[^cond]: 49 coefficients on 384 persistent observations give
`cond(X'X) ≈ 2.2e9`, so the normal equations admit ~`5e-7` relative error in
*any* implementation. Measured difference is `6.5e-8`. See
[`PORTING_LOG.md`](PORTING_LOG.md).

Not yet ported: the `Figure/` plotting layer and the text-table printing
functions. See [Scope](#scope) below.

## Quick start

```julia
using VARToolbox

VAR, VARopt = VARmodel(ENDO, 4)          # 4 lags, constant included

VARopt = @set VARopt.nsteps = 24
VARopt = @set VARopt.ident  = "short"    # recursive identification

IR,  VAR = VARir(VAR, VARopt)            # capture the updated model
VD,  VAR = VARvd(VAR, VARopt)
HD,  VAR = VARhd(VAR, VARopt)

IRinf, IRsup, IRmed, IRbar = VARirband(VAR, VARopt)
```

## Three things that will surprise a MATLAB user

### 1. Options are immutable — update them with `@set`

```julia
VARopt = @set VARopt.nsteps = 24         # NOT VARopt.nsteps = 24
```

MATLAB passes structs by value, and the toolbox depends on it: `SR.m` sets
`VARopt.ident = 'sign'` on its local copy and expects the caller's options
to be untouched. A mutable Julia struct would leak that back and silently
corrupt every band and plot computed afterwards. `@set` is re-exported, so
no extra `using` is needed.

### 2. Identification is late-bound, and each call returns a new model

`VARmodel` performs no identification — `B`, `PSI`, `Fp` come back empty.
`VARir`/`VARvd`/`VARhd` fill them in and return an **updated model** as their
second output. Capture it:

```julia
IR, VAR = VARir(VAR, VARopt)             # not just `IR = VARir(...)`
```

### 3. `IR` and `VD` are indexed differently

```julia
IR[step, variable, shock]                # (nsteps × nvar × nshocks)
VD[step, shock, variable]                # (nsteps × nshocks × nvar)
```

The last two dimensions are transposed relative to each other. This is
upstream behaviour, verified numerically (`sum(VD[h, :, k]) == 100`), and
preserved. `HD.shock` follows `VD`'s convention, not `IR`'s.

## Identification schemes

Set `VARopt.ident` before calling. Support differs by function — this
matrix is upstream's, not a port limitation:

| `ident` | `VARir` | `VARvd` | `VARhd` |
|---|---|---|---|
| `"short"` zero contemporaneous | ✓ | ✓ | ✓ |
| `"long"` zero long-run | ✓ | ✓ | ✓ |
| `"sign"` sign restrictions | ✓ | ✓ | ✓ |
| `"iv"` external instrument | ✓ | error | error |
| `"exog"` exogenous variable | ✓ | error | error |

- **`"iv"` instruments the first column of `ENDO`.** Variable order is
  load-bearing; set `VAR.IV` before calling.
- **`"sign"` requires `B` to be set already** — go through `SR`, don't set
  `ident = "sign"` by hand.
- **Combining IV and sign restrictions** is a two-stage sequence: run
  `VARir` with `ident = "iv"` so `Biv` is populated, then call `SR` with a
  `SIGN` matrix that omits that shock's column.

## Reproducibility

Every stochastic routine takes randomness through an explicit channel:

```julia
VARirband(VAR, VARopt; rng = MersenneTwister(42))
```

and can also consume a pre-drawn pool via `DrawSource`, which is how the
test suite replays MATLAB's own draws to compare bootstrap bands exactly
rather than distributionally. The MATLAB original offers neither.

## Deliberate deviations from the MATLAB source

The port reproduces upstream behaviour — including its defects — so that
numerical comparison stays meaningful. Every such case is documented in the
relevant docstring and in [`PORTING_LOG.md`](PORTING_LOG.md). Four places
deviate deliberately:

1. **`SR`'s draw-index off-by-one.** Upstream returns `IR`/`VD`/`HD`
   computed from a *different* draw than the `B` it returns alongside them.
   Replicated by default; pass `fix_draw_index = true` for the corrected
   behaviour.
2. **`OLSmodel(y, x, 3)`** raises a descriptive error instead of dying on a
   dimension mismatch (upstream adds two columns but increments the
   coefficient count by three).
3. **`LPmodel`** raises a descriptive error when `nlag` is unset, instead of
   failing on a missing struct field (`LPoption.m` never defines it).
4. **`LPmodel`'s `EXOG`/`LAGS` dimension checks** are actually reachable
   here; upstream places them after an `error(...)` call.

Renamed because `const` is reserved in Julia: `VAR.const` → `VAR.constant`,
`HD.const` → `HD.constant`, and the `const_` keyword argument.

## Scope

**Included:** `Utils/`, `Auxiliary/`, `OLSmodel`, `VARmodel`, `ARDLmodel`,
the identification layer, error bands, local projections, and the analytical
part of `Stats/`.

**Not included:**

- **`Figure/`** (20 files). Will land as a package extension
  (`ext/VARToolboxPlotsExt.jl`) so core estimation never depends on a
  plotting stack.
- **The text-table printing layer** — `VARprint`, `OLSprint`, `ARDLprint`,
  `TabPrint`, `mprint`, `roundnum2cell`. Display-only; produces no number
  any downstream calculation uses.
- **`SummStats`** and its `adf`/`ztcrit` dependency.
- **`ExportFig/`** — a Ghostscript wrapper superseded by Julia's native
  vector export.

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite needs **no MATLAB licence**: fixtures are committed alongside
the generators that produced them, with a manifest verified before any test
runs. The manifest carries a canonical hash over the fixture *data* — stable
across regeneration, so generator drift fails loudly — plus a byte hash
recorded for integrity. The byte hash is informational only: MATLAB stamps a
creation timestamp into every `.mat` header, so it changes on any
regeneration even when the data is identical.

Regenerating fixtures does need MATLAB:

```bash
matlab -batch "addpath('test/fixtures/gen'); gen_phase1('test/fixtures')"
julia --project=test test/fixtures/update_manifest.jl
```

Four of the generators use *recording shims* under `test/fixtures/gen/shim/`
— copies of upstream files whose only change is routing their random-number
calls through a wrapper that saves each draw. `make_shim.py` regenerates
them in byte mode, so the diff against upstream stays a two- or three-line
audit trail.

## Attribution

Original MATLAB toolbox © Ambrogio Cesa-Bianchi,
<https://github.com/ambropo/VAR-Toolbox>. This port is a derivative work
under GPL-3.0.
