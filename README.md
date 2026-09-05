# VARToolbox.jl

[![CI](https://github.com/hirothreading/VARToolbox.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hirothreading/VARToolbox.jl/actions/workflows/CI.yml)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE)

A Julia port of **Ambrogio Cesa-Bianchi's VAR Toolbox** (MATLAB, v3.1):
reduced-form VAR estimation, five identification schemes, impulse responses,
forecast error variance decompositions, historical decompositions, bootstrap
error bands, local projections and ARDL models.

Licensed **GPL-3.0**, inherited from the original. See
[`PROVENANCE.md`](PROVENANCE.md) for the upstream commit this derives from.

## Status

The estimation and inference pipeline is complete and validated against the
MATLAB original. **2,685 tests pass**, almost all of them exact
golden-master comparisons at `rtol = 1e-10` against output dumped from
MATLAB R2026a.

The port implements the **VAR Toolbox 3.1 API** and carries the **numerical
bug fixes upstream shipped in VAR Toolbox 4.0** (June–August 2026): residual
degrees of freedom in `OLSmodel`/`VARmodel`/`ARDLmodel`, IV sample
alignment with a gappy instrument, `shut`, lagged exogenous terms in the
historical decomposition, the sign-restriction rotation budget, the
posterior-draw residuals and the draw index in `SR`, the wild bootstrap
with several instruments, and more. Each fix is validated against the
v4.0 MATLAB code itself (`test/fixtures/phase9.mat`). The v4.0 *API
redesign* is not ported — see [Upstream versions](#upstream-versions).

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

The plotting layer ships as a package extension activated by `using Plots`
— see [Plotting](#plotting). Every figure the MATLAB toolbox produces for
the four replications is reproduced in
[`examples/figures/`](examples/figures/), rendered from MATLAB-recorded
draws so the curves are numerically identical, and compared side-by-side
against MATLAB R2026a's own rendering.

Not ported: the text-table printing functions and `SummStats`. See
[Scope](#scope) below.

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

## Plotting

The plotting layer is a package extension (`VARToolboxPlotsExt`): it
activates when Plots.jl is loaded, so core estimation never drags in a
plotting stack.

```julia
using VARToolbox, Plots                  # `using Plots` activates the extension

VARopt = @set VARopt.vnames  = ["Unemployment", "Inflation", "Fed Funds"]
VARopt = @set VARopt.figname = "graphics/SW_"

VARirplot(IRbar, VARopt, IRinf, IRsup)   # one PDF per shock: SW_IR_1.pdf, ...
VARvdplot(VDbar, VARopt)                 # SW_VD.pdf
VARhdplot(HD, VARopt)                    # one PDF per variable: SW_HD_1.pdf, ...
```

The six plot functions (`VARirplot`, `VARvdplot`, `VARhdplot`, `SRirplot`,
`SRvdplot`, `SRhdplot`) read `vnames`, `snames`, `figname`, `FigSize` (in
cm, as in MATLAB), `quality`, `pick`, `firstdate` and `frequency` off
`VARopt`, write `figname`-prefixed PDFs exactly like the MATLAB originals,
and additionally return their `Plots.Plot` objects. The `Figure/` primitives
they stand on are exported too: `PlotSwathe`/`PlotSwathe!`, `AreaPlot!`,
`cmap`, `DatesPlot!`, `FigSize`, `SaveFigure`.

The fidelity target is v3dot1's plot layer with v3dot0's `quality == 2`
save branches restored — v3dot1 dropped them while defaulting
`quality = 2`, so the stock v3dot1 layer silently saves nothing. That and
every other deviation (marker glyph, legend geometry, all quality levels
writing vector PDFs) is listed in [`PORTING_LOG.md`](PORTING_LOG.md)
("Phase 8").

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
  load-bearing; set `VAR.IV` before calling. The instrument may be `NaN`
  where unobserved — leading, trailing or interior gaps are all handled by
  date (v4.0 fix; v3.1 mishandled trailing and interior gaps).
- **`method = "bs"` is rejected under `"iv"`** (as in v4.0): the i.i.d.
  residual bootstrap breaks the instrument/innovation pairing. Use
  `"wild"`.
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

## Upstream versions

Upstream (<https://github.com/ambropo/VAR-Toolbox>) released **v4.0** in
June 2026 and a bug-fix commit in August 2026. v4.0 is an API redesign —
`VARmodel(ENDO, nlag, const, VARopt, …)` performs identification and
inference itself, and `VARir`/`VARvd`/`VARhd`/`VAR*band` no longer exist.
This package keeps the **v3.1 API** it was ported from, but adopts every
v4.0 fix to a computation the two versions share. The full list, with what
was deliberately *not* adopted (the API redesign, the moving block
bootstrap, single-pass `sign+iv`, narrative restrictions, LP-IV, the
`const = 3` removal, an undocumented IV sign convention), is in
[`PORTING_LOG.md`](PORTING_LOG.md), "Phase 9".

Behaviour that changed relative to VAR Toolbox 3.1 as a result:

| function | v3.1 (was replicated) | now (v4.0) |
|---|---|---|
| `OLSmodel` | `tprob`, `bint` on `nobs` df | `nobs - nvar` |
| `VARmodel` | `eq[j].tprob` on `nobs - nvar*nlag` df | `nobs - ntotcoeff` |
| `ARDLmodel` | `bint` on `nobse`, `tprob` on raw `nobs` | `nobse - nvar` |
| `CommonSample` | `lo` counted every non-leading NaN row | trailing run only; empty input OK |
| `VARir`, `ident="iv"` | last-`n`-rows residual alignment | date mask |
| `VARir`, `shut != 0` | Wold recursion ignored `shut` after impact | companion recursion |
| `VARhd`, `nlag_ex > 0` | lagged exogenous terms dropped | included; `HD.endo` reconstructs `Y` |
| `VARirband`, `wild` + `iv` | `k` instrument signs summed into the residual weight | one weight per observation |
| `VARirband`, `bs` + `iv` | ran | rejected |
| `SignRestrictions` | `sr_rot - 1` rotations; `sr_rot = 1` always failed | exactly `sr_rot` |
| `SR`, `sr_mod = 1` | drawn model kept OLS residuals | `resid = Y - X*Ft_draw` |
| `SR` | `IR`/`VD`/`HD` from the draw *after* `B`'s | same draw (`fix_draw_index = false` restores v3.1) |
| `LPmodel`, no `LAGS` | no intercept at all | deterministics always included |
| `MovAvgCent`, `MovCorrCent` | all-NaN on too-short input | error |
| `PairCorr` | unguarded logs | warn and skip on non-positive levels |
| `DatesPlot`, monthly | `0.0833` step, fractional ticks | exact `1/12` |

## Deliberate deviations from the MATLAB source

The port reproduces upstream behaviour — including its defects — so that
numerical comparison stays meaningful. Every such case is documented in the
relevant docstring and in [`PORTING_LOG.md`](PORTING_LOG.md). Beyond the
v4.0 fixes above, three places deviate deliberately:

1. **`OLSmodel(y, x, 3)`** raises a descriptive error instead of dying on a
   dimension mismatch (upstream adds two columns but increments the
   coefficient count by three).
2. **`LPmodel`** raises a descriptive error when `nlag` is unset, instead of
   failing on a missing struct field (`LPoption.m` never defines it).
3. **`LPmodel`'s `EXOG`/`LAGS` dimension checks** are actually reachable
   here; upstream places them after an `error(...)` call.

Renamed because `const` is reserved in Julia: `VAR.const` → `VAR.constant`,
`HD.const` → `HD.constant`, and the `const_` keyword argument.

## Scope

**Included:** `Utils/`, `Auxiliary/`, `OLSmodel`, `VARmodel`, `ARDLmodel`,
the identification layer, error bands, local projections, the analytical
part of `Stats/`, and the plotting layer (the six `VAR/*plot` functions
plus the `Figure/` primitives they use — `PlotSwathe`, `AreaPlot`, `cmap`,
`DatesPlot`, `FigSize`, `SaveFigure` — as the `VARToolboxPlotsExt`
extension).

**Not included:**

- **The rest of `Figure/`** — `BarPlot`, `PlotLine` (broken upstream: it
  calls a function that does not exist in the toolbox), `PlotStatesShaded`,
  `pantone`, `LegPlot`, and `cmap`'s `parula` tail for more than 9 series.
  None is used by any toolbox figure.
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

The test suite needs **no MATLAB licence** and no files outside this
package: the replication `.xlsx` sheets are vendored byte-for-byte in
`data/` (checked against the MATLAB originals whenever the monorepo is
present), and fixtures are committed alongside the generators that produced
them, with a manifest verified before any test runs. The manifest carries a canonical hash over the fixture *data* — stable
across regeneration, so generator drift fails loudly — plus a byte hash
recorded for integrity. The byte hash is informational only: MATLAB stamps a
creation timestamp into every `.mat` header, so it changes on any
regeneration even when the data is identical.

Regenerating fixtures does need MATLAB and a checkout of the fork:
`gen_phase1..7.m` read its `v3dot1/`, `gen_phase9.m` (the v4.0 fixes) reads
its `v4dot0/`.

```bash
matlab -batch "addpath('test/fixtures/gen'); gen_phase1('test/fixtures')"
julia --project=test test/fixtures/update_manifest.jl
```

Four of the generators use *recording shims* under `test/fixtures/gen/shim/`
— copies of upstream files whose only change is routing their random-number
calls through a wrapper that saves each draw. `make_shim.py` regenerates
them in byte mode, so the diff against upstream stays a two- or three-line
audit trail.

## Attribution and citing

Original MATLAB toolbox © Ambrogio Cesa-Bianchi,
<https://github.com/ambropo/VAR-Toolbox>. This port is a derivative work
under GPL-3.0. See [`PROVENANCE.md`](PROVENANCE.md) for the exact upstream
commit and [`PORT_SUMMARY.md`](PORT_SUMMARY.md) for how the port was
validated.

If you use this package in academic work, please cite the original
toolbox (Cesa-Bianchi, A., *VAR Toolbox*,
<https://github.com/ambropo/VAR-Toolbox>) alongside this port
(<https://github.com/hirothreading/VARToolbox.jl>).
