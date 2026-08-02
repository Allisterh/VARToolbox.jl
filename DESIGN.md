# VARToolbox.jl — design decisions

Companion to `../JULIA_PORT_PROMPT.md` (revision 2). That brief settles
three decisions (D1–D3) and sets three technical mandates (M1–M3). This
document records how they are actually implemented, plus the module layout
and fixture mechanics.

---

## D1 — Naming: MATLAB-verbatim functions, suffixed types

Exported **functions keep their MATLAB names exactly**: `VARmodel`,
`VARoption`, `VARir`, `VARvd`, `VARhd`, `VARirband`, `VARvdband`, `SR`,
`SignRestrictions`, `VARdrawpost`, `OLSmodel`, `ARDLmodel`, `LPmodel`,
`LPoption`, `LPmakex`, `CommonSample`, `VARmakexy`, `VARmakelags`, `getqr`,
`OrthNorm`, `trimr`, `L`, `seqa`, `ptrend`, `tdis_prb`, `tdis_inv`,
`fdis_prb`, `beta_pdf`, `beta_inv`, and the `nan*` family.

**Types take distinct suffixed names**: `VARModel`, `VAROptions`, `LPModel`,
`LPOptions`, `ARDLModel`, `OLSOutput`, `SROutput`, `EquationResult`.

The reason is concrete rather than stylistic. The MATLAB idiom is

```matlab
[VAR, VARopt] = VARmodel(ENDO, nlag);
```

If the Julia type were named `VAR`, the binding `VAR = VARmodel(...)`
shadows it for the rest of the scope and the constructor becomes
unreachable. Suffixing costs nothing in recognisability.

### Names deliberately *not* exported

`cols`, `rows`, `columns`, `vec`, `lag`. These are MATLAB shims for things
Julia already has (`size(X,2)`, `size(X,1)`, `Base.vec`) or are too generic
to claim in a user's namespace. They exist as unexported internals
(`_cols`, `_rows`, `_vec`, `lag`) so call sites stay visually close to the
MATLAB source. Reach them as `VARToolbox.lag` etc.

### Forced rename

MATLAB's `OLS.const` / `VAR.const` becomes **`constant`**. `const` is
reserved in Julia and cannot be a struct field name. This is the only rename
imposed by the language rather than by choice.

---

## D2 — Fidelity: replicate, document, and pin

Where the MATLAB source is odd, the oddity is reproduced, given a docstring
admonition, and pinned with a test. Silent improvement is not permitted;
divergence is a later, separate, deliberate PR.

Two carve-outs, both recorded in `PORTING_LOG.md`:

1. **`VARmodel`'s `const` docstring bug** — documentation says default `0`,
   code defaults to `1`. The *documentation* is fixed; behaviour stays `1`.
   (Agreed in the brief.)
2. **`OLSmodel`'s `const == 3`** — the upstream branch is internally
   inconsistent (adds two columns, increments `nvar` by three) and always
   raises a dimension-mismatch error. Replicating a guaranteed crash helps
   nobody, so the port raises an explicit `ArgumentError` naming the defect.
   The fixture records that MATLAB errors here, so the behavioural contract
   ("this input is rejected") is still validated.

Everything else found so far is replicated verbatim: `CommonSample`'s `lo`
miscount, `nanmedian`'s zero-median-to-NaN conversion, `tdis_prb`'s clamp,
`tdis_inv`'s precision collapse near `p = 0.5`, `L`'s inverted sign
convention, `OLSmodel`'s `nobs` degrees of freedom. See `PORTING_LOG.md`.

---

## D3 — Options structs are immutable

`VAROptions` and `LPOptions` are **immutable**, updated by reconstruction
via Setfield's `@set` (re-exported, so no extra `using` is needed):

```julia
VARopt = VARoption()
VARopt = @set VARopt.nsteps = 24
VARopt = @set VARopt.ident  = "short"
```

This is not stylistic preference. MATLAB passes structs by value and the
toolbox depends on it:

- `SR.m:94` and `SR.m:178` set `VARopt.ident = 'sign'`
- `SignRestrictions.m:100` sets `VARopt.nsteps = VARopt.sr_hor`

Under MATLAB's copy semantics the caller's `VARopt` is untouched. With a
`mutable struct` the caller would get back `ident = "sign"` and
`nsteps = 1`, silently corrupting every band and plot computed afterwards.
Immutability reproduces the MATLAB contract exactly and makes the failure
mode impossible.

The cost is that `@set` is unfamiliar to MATLAB users. It is documented
prominently in the README as the single biggest surface difference.

---

## M1 — `prctile` is Hazen (type 5), not Julia's default

Every MATLAB `prctile` call routes through `VARToolbox._prctile`, which is
`quantile(v, p/100; alpha = 0.5, beta = 0.5)` with NaN omission.

Julia's default `quantile` is type 7 and differs by percent-level amounts at
`pctg = 95` — the toolbox's own default confidence level. A test asserts
that the naive translation does **not** match, so this cannot regress into a
silently-wrong error band.

Call sites: `SR.m:168,172`, `VARirband.m:177-179`, `VARvdband.m:162-164`,
`Stats/censor.m:30-31`, `Auxiliary/winsor.m:35`.

---

## M2 — Injectable randomness

Every stochastic routine takes randomness through an explicit channel:

```julia
f(...; rng::AbstractRNG = Random.default_rng(), draws = nothing)
```

A non-`nothing` `draws` is consumed instead of sampling. This exists so that
MATLAB's own draws can be dumped to a fixture and fed into the Julia
function, turning bootstrap and sign-restriction validation into **exact**
golden-master checks at `rtol = 1e-10` rather than distributional
comparisons at percent-level tolerance.

The full stochastic surface is eight call sites in five files:
`VARirband.m:84,88,92`, `VARvdband.m:80,83`, `SignRestrictions.m:64,92`,
`OrthNorm.m:19`, `VARdrawpost.m:43,53`.

Side benefit: the Julia package becomes reproducible and safely
parallelisable, which the MATLAB original is not.

---

## M3 — Hand-rolled distributions are replicated, not substituted

`Auxiliary/`'s distribution functions are built on base MATLAB's `betainc`
(`= cdf(Beta(a,b), x)`) and carry **no Statistics Toolbox dependency**.
Their numerical quirks exceed the golden-master tolerance, so they are
replicated rather than replaced with `Distributions.jl` closed forms. Detail
in `PORTING_LOG.md` §Phase 1.

---

## Module layout

A single flat namespace, files `include`d by `src/VARToolbox.jl`:

```
src/VARToolbox.jl      module, exports, design notes
src/auxiliary.jl       Auxiliary/  — _prctile, distributions, nan* family
src/utils.jl           Utils/      + VARmakexy, VARmakelags, L
src/ols.jl             OLSmodel
src/var.jl             VARmodel, VARoption, VARModel, EquationResult   [Phase 2]
src/ardl.jl            ARDLmodel                                        [Phase 2]
src/stats.jl           Stats/ (analytical part; display layer not ported)
src/identification.jl  VARir, VARvd, VARhd, SR, SignRestrictions, VARdrawpost [Phase 3]
src/bands.jl           VARirband, VARvdband                             [Phase 4]
src/lp.jl              LPmodel, LPoption, LPmakex                       [Phase 5]
src/plotstubs.jl       zero-method stubs + docstrings for the plot layer [Phase 8]
```

**No submodules per folder.** The MATLAB API is one flat namespace;
submodules would force qualified names on readers transferring from the
manual. The folder boundaries survive as file boundaries.

**Plotting is a package extension** (Phase 8):

```
ext/VARToolboxPlotsExt/VARToolboxPlotsExt.jl   entry: imports the stubs, includes the two files
ext/VARToolboxPlotsExt/primitives.jl           Figure/  — cmap, PlotSwathe, AreaPlot, DatesPlot,
                                               FigSize, SaveFigure, legend-row helper
ext/VARToolboxPlotsExt/varplots.jl             VAR/     — VARirplot, VARvdplot, VARhdplot, SR*plot
```

`Project.toml` declares Plots under `[weakdeps]`/`[extensions]`, so the
extension activates on `using Plots` and core estimation never drags in a
plotting stack; no estimation function calls into plotting. The function
*names* live in `src/plotstubs.jl` as zero-method stubs carrying the
docstrings, so `?VARirplot` works (and points at `using Plots`) even
before the extension loads; the extension `import`s and adds the methods.
Plots also sits in `test/Project.toml`, which is what loads the extension
under `Pkg.test` (verified on 1.10 and 1.12). The backend is GR, headless
via `ENV["GKSwstype"] = "100"`; `FigSize` converts MATLAB's centimetre
figure sizes at 1 cm = 72/2.54 px, which GR maps 1:1 to PDF points.

### Replacing the `eval`-built fields

The MATLAB source has 17 `eval` sites in two distinct patterns:

- **16 per-equation** — `VARmodel.m:107-121` (writes `VAR.eq1.beta`, …) and
  `VARprint.m:129,140,151,175,176` (reads them). → `VARModel.eq` as a
  `Vector{EquationResult}`, indexed `VAR.eq[j]`.
- **1 per-horizon** — `LPmodel.m:132` builds `h1`, `h2`, … labels. A
  different pattern, given its own container in Phase 5.

---

## Fixture mechanics

```
test/fixtures/gen/gen_phaseN.m   committed MATLAB generators (functions, not scripts)
test/fixtures/phaseN.mat         committed output, -v7 format
test/fixtures/MANIFEST.toml      two hashes per fixture + generator + MATLAB version
test/fixtures/content_hash.jl    canonical hash over fixture DATA
test/fixtures/update_manifest.jl regenerates the manifest
```

Both generators and their output are committed, so `] test` runs on a
machine with **no MATLAB licence**.

The manifest carries two hashes, because they answer different questions:

- **`content_sha256`** — a canonical hash over the data, stable across
  regeneration. This is what the test suite asserts, and what makes a
  generator that has drifted from its committed output fail loudly.
- **`sha256`** — a hash of the file bytes. Recorded for integrity only, and
  checked as an informational message rather than an assertion, because
  MATLAB stamps a creation timestamp into every `.mat` header. Regenerating
  `phase1.mat` from an unchanged generator changes this hash while leaving
  all 115 fields bit-identical, so asserting on it would fail every honest
  regeneration while still catching nothing that matters.

Generators are **functions**, not scripts, because the `Replic/*/GO_*.m`
scripts open with `clear all`, which wipes the workspace of any script-level
driver that loops over them.

Regeneration:

```bash
matlab -batch "addpath('test/fixtures/gen'); gen_phase1('test/fixtures')"
julia --project=test test/fixtures/update_manifest.jl
```

## Validation policy

| class | method | tolerance |
|---|---|---|
| deterministic | golden-master `.mat`, elementwise | `rtol = 1e-10` |
| stochastic with injected draws | golden-master with MATLAB's draws fed in | `rtol = 1e-10` |
| stochastic without injection | mean/std/quantiles over ≥10,000 draws | ~1–2%, **labelled supplementary** |

Any tolerance looser than `1e-8` carries a written justification at the call
site. Comparisons use `approxeq`, which treats `NaN` as matching `NaN` —
necessary because NaN is meaningful data here (`L` pads with it, `HD`
NaN-pads its first `nlag` rows, the `nan*` family returns it).
