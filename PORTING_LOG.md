# Porting log

One entry per phase. Each records what was validated against MATLAB, the
tolerance used, deliberate deviations, and — importantly — **what was not
checked**. Gaps are stated, not skipped.

Golden-master tolerance is `rtol = 1e-10` throughout unless noted.

---

## Phase 0 — blockers and harness

**Validated:** that the MATLAB source runs at all under the reference
environment, since every fixture depends on it.

| script | result |
|---|---|
| `Primer/VARToolbox_Primer.m` | **runs** after one-line path patch; 31 figures written |
| `Replic/SW2001/GO_SW2001.m` | **runs** clean (31 s) |
| `Replic/BQ1989/GO_BQ1989.m` | numerics complete; **fails at final `SaveFigure`** |
| `Replic/Uhlig2005/GO_Uhlig2005.m` | numerics complete; **fails at final `SaveFigure`** |
| `Replic/GK2015/GO_GK2015.m` | numerics complete; **fails at final `SaveFigure`** |

**Findings:**

1. `Primer/VARToolbox_Primer.m:6` hardcodes `/Users/jambro/GDrive/VAR-Toolbox/v3dot1`.
   Patched in a scratch copy — `v3dot1/` is never written to. This was the
   one true blocker.
2. The three `SaveFigure` failures are `SaveFigure(...,1)`, the `quality=1`
   Ghostscript path, which calls `export_fig` — bundled in `v3dot0/ExportFig`
   but **dropped from `v3dot1`**. This is the documented gotcha, confirmed.
   All three are the *last statement* of their script, so every numerical
   result is produced before the failure. Not a blocker for fixtures, since
   generators do no plotting.
3. The `Replic/*/GO_*.m` scripts open with `clear all`, which wipes the
   workspace of any driver that loops over them. Fixture generators are
   therefore written as **functions**, not scripts.
4. `Replic/*/GO_*.m` do **not** hardcode paths (contrary to the brief's
   assumption) — they rely on the toolbox already being on the path.

**Harness:** `test/fixtures/gen/*.m` (committed generators) →
`test/fixtures/*.mat` (committed output) → `MANIFEST.toml` (SHA-256 per
fixture). `runtests.jl` verifies hashes before running anything, so a
fixture that drifted from its generator fails loudly. The suite needs no
MATLAB licence.

**Not checked:** the `quality=1` export path (out of scope with plotting);
whether the Primer's *numbers* are correct, only that it completes.

---

## Phase 1 — Utils + Auxiliary (+ the OLS engine)

**Validated against:** `CommonSample`, `VARmakelags`, `VARmakexy` (all four
`const` values), `lag`, `trimr`, `L`, `seqa`, `ptrend`, `getqr`, `Num2NaN`,
`NaN2Num`, `prctile`, `tdis_prb`, `tdis_inv`, `beta_pdf`, `beta_inv`,
`fdis_prb`, `nanmean`, `nansum`, `nanstd`, `nanvar`, `nansem`, `nanmedian`,
`OLSmodel` (`const` 0/1/2, all 17 output fields).

**Tolerance:** `rtol = 1e-10`. Achieved agreement is substantially tighter
— `tdis_inv` matches to 1.7e-15 relative.

**Result:** 125 / 125 tests pass.

### M1 — `prctile` (confirmed and handled)

MATLAB's `prctile` is the Hazen / type-5 definition; Julia's default
`quantile` is type 7. Measured on the fixture vector:

| p | MATLAB | Julia default | Julia `alpha=beta=0.5` |
|---|---|---|---|
| 2.5 | 1.4000 | 1.54 | 1.4000 |
| 95 | 9.6000 | 9.32 | 9.6000 |

All call sites route through `VARToolbox._prctile`. A test asserts that the
*naive* translation does **not** match, so the distinction cannot silently
regress.

### M3 — hand-rolled distributions (confirmed; replicated, not substituted)

`tdis_prb`/`tdis_inv`/`beta_inv`/`beta_pdf`/`fdis_prb` use base MATLAB's
`betainc` — **no Statistics Toolbox dependency**, contrary to what the
Auxiliary headers imply. Two of them could not be replaced by
`Distributions.jl` closed forms without breaking agreement:

- **`tdis_prb` clamps** `x2` to `1 - 1e-12`. Since `Beta(n/2, 1/2)` is
  singular at 1, this perturbs the result by ~1e-6 near `x = 0` — four
  orders of magnitude above the test tolerance. Clamp replicated.
- **`tdis_inv` loses ~7 digits near `p = 0.5`.** `beta_inv` stops at a
  Newton residual of ~1e-14 and the final `sqrt` amplifies it. MATLAB
  returns `3.2148354647791181e-07` where the exact answer is `0`; this port
  reproduces that to 16 significant digits. Harmless for `OLSmodel`'s
  `tdis_inv(0.025, nobs)`, but pinned by test.

### Upstream defects found, replicated, and pinned by test

1. **`CommonSample`'s `lo` is not the trailing-NaN count.** Its loop never
   stops at the first non-NaN, so `lo` counts *every* NaN row outside the
   leading run. Coincides with the documented meaning only when the data has
   no interior gaps. The fixture deliberately includes an interior NaN row;
   the test asserts `lo == 3`, not `2`.
2. **`nanmedian` turns an exactly-zero median into `NaN`** (`i = find(y==0);
   y(i) = i + nan;`). Replicated and pinned.
3. **`OLSmodel`'s `const == 3` is broken.** It prepends two columns
   (`[ones trend.^2]`) but increments `nvar` by three, so `xpxi` is built at
   the wrong dimension and MATLAB errors. Note also that `VARmakexy`'s
   `const == 3` prepends *three* columns — the two functions disagree about
   what `3` means. **Deviation:** rather than replicate a guaranteed crash,
   `OLSmodel` throws an `ArgumentError` naming the defect. The fixture
   records that MATLAB errors here.
4. **`OLSmodel`'s `tprob` uses `nobs` degrees of freedom, not `nobs - nvar`.**
   Replicated. `VARmodel` recomputes its own `tprob` with the correct
   `nobs - ncoeff`, so per-equation `tprob` inside a fitted VAR differs from
   `OLSmodel`'s — that is upstream behaviour, not a port artefact.
5. **`L`'s sign convention is inverted from its name**: positive `nlag` is a
   *lead*. Preserved, with a docstring warning.
6. **`Utils/columns.m` declares `function c = cols(DATA)`** — filename and
   function name disagree. MATLAB dispatches on the filename so it works.

### Deliberate deviations

- **`cols`, `rows`, `columns`, `vec` are not exported.** They are MATLAB
  shims for `size(X,2)`, `size(X,1)` and `Base.vec`. Kept as unexported
  internals (`_cols`, `_rows`, `_vec`) so the package does not claim such
  generic names in a user's namespace.
- **`lag` is not exported** for the same reason; reachable as
  `VARToolbox.lag`.
- **`OLSOutput.const` renamed to `constant`** — `const` is reserved in Julia.
  This rename propagates to `VARModel` in Phase 2.
- **`nanmovavg` and `nancumsum` are not ported.** `nancumsum` is referenced
  nowhere outside `Auxiliary/`. `nanmovavg` is broken upstream: its inner
  loop reads from the output array `y` (all `NaN` at that point) rather than
  the input `x`, so it returns all-`NaN` for every input, and its `aux_x`
  local is assigned but never used. Porting a function that cannot return a
  non-NaN value serves no one.
- **`nanmax`/`nanmin`** are legacy shims that forward to `max`/`min`;
  provided for API parity only.

**Not checked:** `adf`, `ztcrit`, `winsor`, `censor`, `roundn`, `tdiff`,
`XoX`, `datatreat`, `cell2num`, `roundnum2cell`, `NoSpace`, `TabPrint`,
`mprint`, `m2tex`, `ws2struct`, `clearex`, `SupTitle`, `rgb`, `hatchfill2`.
Of these, only `roundnum2cell`/`TabPrint`/`mprint` are reachable from the
core call graph, and only from the *printing* functions (`VARprint`,
`OLSprint`, `ARDLprint`) — deferred to Phase 2. The rest are either
plotting-layer (deferred with `Figure/`) or unreferenced by the core.

---

## Phase 2 — core estimation

**Validated against:** `VARoption` (all 20 defaults), `VARmodel` (6
lag/deterministic configurations plus all 3 `nlag` vs `nlag_ex` branches
plus real Stock-Watson data), `ARDLmodel` (3 configurations).

**Tolerance:** `rtol = 1e-10`. **Result:** 581 / 581 tests pass.

Fixtures use real SW2001 data at exactly the configuration `GO_SW2001.m`
uses (`nlag = 4`, `const = 1`), so a Phase 6 end-to-end discrepancy can be
localised to identification rather than estimation. Agreement there is 2e-12
absolute on a coefficient matrix spanning 0.0042 to 1.97.

### The 17 `eval` sites

16 of them are gone: `VARModel.eq::Vector{EquationResult}`, indexed
`VAR.eq[j]`, replaces the 11 `eval` writes at `VARmodel.m:107-121` and the 5
reads at `VARprint.m:129-176`. The 17th (`LPmodel.m:132`, per-horizon
labels) is a different pattern and is handled in Phase 5.

### Upstream behaviour replicated and documented

1. **`VARmodel` recomputes `tprob` on `nobse - ncoeff` dof**, where `ncoeff
   = nvar*nlag` counts only lag coefficients — excluding the constant and
   any exogenous regressors. `OLSmodel` uses `nobs`. The two disagree by
   construction; the `VARmodel` value is what `VAR.eq[j].tprob` carries.
2. **`ARDLmodel` uses two different dof conventions internally**:
   `tdis_inv(.025, nobse)` for `bint` but `tdis_prb(tstat, nobs)` for
   `tprob`. So its confidence intervals and p-values are not mutually
   consistent. Replicated.
3. **`ARDLmodel.nvar` counts coefficients, not variables** (`ncoeff +
   ncoeff_ex + const`). Misleading name kept for API parity.
4. **`VARmodel`'s `const` docstring bug fixed, behaviour unchanged** — the
   agreed carve-out. Default is 1.

### Notes on numerical fidelity

- `Ft` is computed with `cholesky(Symmetric(X'X))`, matching MATLAB's
  dispatch on symmetric positive-definite systems. Julia's generic `\` uses
  LU with partial pivoting; both are correct but they round differently, and
  matching the algorithm keeps agreement near machine precision.
- **Fixture input data is generated by an explicit LCG, not `rand`.** The
  first draft used smooth trigonometric series, which made each variable
  nearly a linear combination of its own lags and produced designs with
  `RCOND ~ 1e-19`. Comparing two implementations at `rtol = 1e-10` on a
  singular system is meaningless. White-noise-like input keeps the design
  well conditioned, and an explicit LCG keeps the fixture byte-reproducible
  without depending on MATLAB's RNG.
- **Test tolerance policy refined.** `approxeq` now passes an element when
  it is close relatively *or* relative to the array's own scale
  (`atol = rtol * maxabs`). A VAR coefficient matrix spans four orders of
  magnitude; demanding that the smallest coefficient agree to `rtol` in its
  own units is a test of the linear-algebra backend, not of the port. The
  scaled-atol form still catches a single element that is genuinely wrong,
  which a norm-based `isapprox` would hide.

### Deferred within Phase 2

**`Stats/` and the print functions are not yet ported.** `VARprint`,
`OLSprint`, `ARDLprint` (and their dependencies `TabPrint`, `mprint`,
`roundnum2cell`) are a display layer with no effect on any numerical result;
`Stats/` (`SummStats`, `CorrTable`, `MovAvg`, `MovCorr`, `PairCorr`,
`CrossCorrelation`, …) is independent of the VAR pipeline entirely. Both are
scheduled after Phase 5 so that identification, error bands and local
projections — the parts users actually depend on — land first.

Note for whoever picks up `Stats/`: `SummStats.m` calls `skewness` and
`kurtosis`, which **are** Statistics Toolbox functions. That is an
additional toolbox dependency beyond the six `prctile`/`norminv` sites
catalogued in the brief.

---

## Phase 3 — identification and IR/VD/HD

**Validated against:** `VARir` (all five `ident` schemes, plus `recurs`,
`impact` and `shut` variants), `VARvd`, `VARhd` (including the exogenous
branch), `SignRestrictions` (`sr_hor = 1` and `sr_hor = 4`), `SR` end to
end, `VARdrawpost`.

**Tolerance:** `rtol = 1e-10`, except the Gertler-Karadi IV block at `1e-6`
(justified below). **Result:** 740 / 740 tests pass.

### M2 realised — sign restrictions validated exactly, not distributionally

The brief's plan was to compare bootstrap and rotation output
distributionally at 1–2%. Instead, `gen/shim/SignRestrictions.m` — a copy of
the upstream file whose **only** change is routing its two `randn(` call
sites through a recording wrapper (`diff` confirms two lines) — saves the
exact draw sequence MATLAB consumed. The Julia side replays it through a
`DrawSource`.

Consequence: `SignRestrictions` and the whole of `SR` (with `sr_mod = 0`,
which removes the only other stochastic step) are **exact golden-master
checks at `rtol = 1e-10`**, covering `Ball`, `IRall`, `VDall`, the median
and percentile bands, and the historical decompositions. A test also asserts
the Julia run consumed exactly as many draws as MATLAB did, which confirms
both implementations accepted and rejected the same rotations.

`gen/shim/make_shim.py` regenerates the shim in byte mode, preserving the
upstream CRLF line endings so the two-line diff stays auditable.

### Layout finding: IR and VD are indexed differently

`CLAUDE.md` states that `IR` is `(nsteps × nvar × nshocks)` and that "`VD`
matches". **It does not.** Verified numerically:

- `IR[step, variable, shock]` — `IR[1, :, 1]` equals column 1 of `B`.
- `VD[step, shock, variable]` — `sum(VD[h, :, k]) == 100` for every horizon
  and variable, while summing over the third index gives 120.9 and 46.7.

`HD.shock` follows `VD`'s convention, not `IR`'s. All three are pinned by
test and documented with admonitions on the respective functions.

### Upstream defect: `SR` returns `IR`/`VD`/`HD` from the wrong draw

`SR.m` files each draw's model under a label built from `jj` **before** `jj`
is incremented, so accepted draw `k` is stored as `"draw(k-1)"`. At the end
it selects `sel` — the accepted draw whose `B` is closest to the median —
and reads back `VAR_draw.("draw" + sel)`, which is draw `sel+1`.

Verified on Stock-Watson data with `sr_mod = 0` (so all draws share one
model and the mapping is unambiguous): **`SRout.B` was accepted draw 24
while `SRout.IR` matched accepted draw 25 to 0.0e+00.** So the impulse
responses, variance decomposition and historical decomposition that `SR`
returns as its point estimate do not correspond to the `B` it returns
alongside them.

There is also a latent crash: when `sel == ndraws`, the label
`"draw(ndraws)"` was never created and MATLAB errors on a non-existent
field. Probability roughly `1/ndraws` per call.

Replicated by default so golden-master comparison stays meaningful, with
`fix_draw_index = true` available to use the draw that actually corresponds
to `B` (and which also removes the crash). Both behaviours are pinned by
test. The `med`/`inf`/`sup` bands are computed across all draws and are
unaffected.

### Bug found in this port by the fixtures

`_wold` initially allocated `Fp` with `nlag` slices, mirroring MATLAB's
`zeros(nvar,nvar,nlag)` preallocation. But MATLAB's else-branch assigns
`Fp(:,:,ii) = zeros(nvar,nvar)` for `ii > nlag`, which **grows** the array
to `nsteps`, and the Wold recursion indexes up to `nsteps-1`. Allocating
`max(nlag, nsteps)` fixed it. Caught immediately by the golden-master
comparison of `VAR.Fp`.

### Tolerance justification — the IV block

The Gertler-Karadi specification is 12 lags of 4 variables plus a constant:
49 coefficients on 384 observations of persistent macro series. `cond(X)` is
about 4.7e4, so `cond(X'X)` is about 2.2e9 and solving through the normal
equations admits relative error of order `cond(X'X)·eps ≈ 5e-7` in any
implementation.

Measured difference against MATLAB is **6.5e-8 norm-relative** (3.1e-6
absolute on coefficients up to 32.1) — an order of magnitude inside that
bound. Both results are equally valid solutions of the same ill-conditioned
system, so this block runs at `rtol = 1e-6`. The test asserts
`cond(X'X) > 1e8` alongside, so the justification is self-checking: if the
design ever becomes well conditioned the loosened tolerance stops being
warranted.

Note this conditioning is inherent to what the toolbox does, not to the
fixture: `GO_GK2015.m` uses exactly this specification.

### Other findings

- **`GO_GK2015.m` applies `Num2NaN` before selecting columns.** The sheet
  encodes missing values as the sentinel `123456789` and the `ff4_tc`
  instrument is only populated over part of the sample. A first draft of the
  generator skipped this and produced an IV first stage with
  `RCOND ~ 1.5e-16`. Worth flagging for anyone writing new example scripts.
- **`VARhd`'s exogenous contribution ignores lags of the exogenous
  variables.** It uses only `F[:, nvar*nlag+const+ii]`, the contemporaneous
  coefficient. With `nlag_ex > 0` the lagged exogenous coefficients are
  silently omitted, so `HD.endo` does not reconstruct the data exactly.
  Replicated and documented.
- **`SignRestrictions` returns `nothing` on rotation-budget exhaustion and
  does not raise** — its error path is commented out upstream. Pinned by a
  test using deliberately impossible restrictions. `SR`'s separate
  `sr_draw` budget is the layer that does raise.

**Not checked:** the `whereToStart > 1` branch of `SignRestrictions` — the
path that combines an instrument-identified first column (`VAR.Biv`) with
sign restrictions on the remaining shocks, and the only consumer of the
`randn(dy,1)` draw site. The brief calls this combination out explicitly, so
it is a real gap. Everything it depends on (`ident = "iv"`, `sigma_b`,
`Biv`, and the square-`SIGN` path through the same loop) is validated;
what is untested is the orthogonal-completion step itself.

---

## Phase 4 — error bands

**Validated against:** `VARirband` (standard and wild bootstrap; `short` and
`long` identification; `pctg = 95` and `68`; with an exogenous variable;
and the `iv` + `wild` branch) and `VARvdband` (standard and wild).

**Tolerance:** `rtol = 1e-10`, except the Gertler-Karadi IV block at `1e-6`
for the conditioning reason given in Phase 3. **Result:** 830 / 830 tests
pass.

### M2, second application: bands are exact, not distributional

The brief proposed validating bootstrap bands by comparing summary
statistics over ~10,000 draws at 1–2% relative tolerance. That tolerance is
*wider than the `prctile` defect M1 identifies* — it would have passed a
build in which every band was computed from the wrong quantile definition.

Instead, `gen/shim/VARirband.m` and `gen/shim/VARvdband.m` route their
`rand(` call sites through a recording wrapper (3 and 2 sites respectively;
`diff` confirms only those lines change) and the recorded uniforms are
replayed in Julia through a `DrawSource`. Every band is therefore compared
elementwise at `rtol = 1e-10`, over `INF`, `SUP`, `MED` and `BAR`.

`DrawSource` tracks normal and uniform pools separately, so a pool recorded
from `randn` can never be consumed by a `rand` call site. Exhausting an
injected pool raises rather than falling back to sampling — a golden-master
test that quietly stopped being one would be worse than no test.

### Behaviour preserved

- **`VARirband` has an `iv` wild-bootstrap branch; `VARvdband` does not.**
  With `method = "wild"`, `VARirband` checks `ident == "iv"` and resamples
  the instrument alongside the residuals
  (`Z = [IV(1:nlag,:); IV(nlag+1:end,:).*rr]`); `VARvdband` always takes the
  plain `rand(nobs,1)` path. Replicated. It rarely bites in practice
  because `VARvd` refuses `ident = "iv"` outright, but it does mean the two
  functions consume different numbers of draws for the same options — which
  matters precisely when injecting recorded draws.
- **Draws with `maxEig >= 0.9999` are discarded and redrawn**, so `ndraws`
  counts accepted replications rather than attempts, and there is no cap on
  attempts. An explosive point estimate loops indefinitely. Upstream.
- **`VARvdband.m` hardcodes its progress interval as `10*ww`** instead of
  reading `VARopt.mult` as `VARirband.m` does. Not reproduced — this port
  prints no progress at all — and noted only because it is a real
  inconsistency between two otherwise near-identical functions.
- Bands use `_prctile` (M1), so they inherit MATLAB's Hazen/type-5
  definition rather than Julia's default.

**Not checked:** the rejection path itself (no fixture drives a draw to
`maxEig >= 0.9999`, so the resample-and-retry branch is exercised only
incidentally); `VARvdband` with an exogenous variable.

---

## Phase 5 — local projections

**Validated against:** `LPoption`, `LPmakex`, `LPmodel` (with lags, unit
shock, without lags, with exogenous regressors, and at `pctg = 68`).

**Tolerance:** `rtol = 1e-10`. No RNG is involved. **Result:** 1050 / 1050
tests pass.

Kept as the independent API it is upstream: its own options struct, its own
result type, no interaction with `VARir`/`VARvd`/`VARhd`.

### The 17th `eval` site

`LPmodel.m:132` builds `LPout.h1`, `LPout.h2`, … one horizon at a time. This
is the only `eval` site that is *not* the per-equation pattern, and it is
now `LPOutput.h::Vector{HorizonResult}`, indexed `LPout.h[hh]`. **All 17
`eval` sites in the toolbox are now eliminated.**

### Upstream defects found

1. **`LPoption.m` does not define `nlag` or `const`, but `LPmodel.m` reads
   both on its 8th and 9th lines.** Calling `LPmodel` with an unmodified
   options struct fails in MATLAB with a bare "Reference to non-existent
   field 'nlag'". The file's own guards (`if ~exist('nlag','var')`,
   `if ~exist('const','var')`) are dead code — the variables were assigned
   from `LPopt` two lines earlier, so they can never fire.

   **Deviation:** this port declares both fields. `nlag` defaults to
   `nothing` and raises an `ArgumentError` naming the problem, preserving
   the contract that it is mandatory while replacing an obscure
   field-access failure with an actionable message. Pinned by a test that
   also asserts the MATLAB struct genuinely lacks the fields.

2. **`bstd_HW` is assigned twice and holds Newey-West values.**

   ```matlab
   LPout.(aux{1}).bstd_HW = OLSout.bstd_HW;
   LPout.(aux{1}).bstd_HW = OLSout.bstd_NW;   % overwrites the previous line
   ```

   So the field named for White standard errors carries Newey-West ones,
   and the White values are lost. Replicated — `HorizonResult.bstd_HW`
   matches MATLAB — with the discarded value also exposed as
   `bstd_White`. Tests pin that the stored field equals `bstd_NW`, that
   `bstd_White` equals `bstd_HW`, and that the two genuinely differ.

3. **The bands do not use the field the model stores.** `INF`/`SUP` are
   computed from `OLSout.bstd_HW` — the *White* errors — while the stored
   `bstd_HW` field holds Newey-West. So the reported bands cannot be
   reconstructed from the reported standard errors. Pinned by test.

4. **`EXOG` and `LAGS` row-count checks are unreachable.** Both are nested
   inside `if nobsS ~= nobsENDO` and placed *after* its `error(...)` call.
   **Deviation:** this port checks all three dimensions properly. This
   cannot change the result of any input MATLAB accepts — it only rejects
   inputs MATLAB would have silently mis-handled.

### Behaviour worth knowing

- **With `LAGS` empty the projection has no intercept at all**, regardless
  of `const_`. The regression is run as `OLSmodel(Y, X, 0)` and the constant
  enters only through `VARmakexy(LAGS, nlag, const_)`. Pinned by a test
  showing `const_ = 0` and `const_ = 1` give identical results in that case.
- **Bands are asymptotic Normal, not bootstrapped**, and carry no correction
  for the serial correlation that overlapping horizons induce by
  construction. They are therefore anti-conservative at longer horizons.
  Upstream's choice; documented on `LPmodel` so a user does not mistake
  them for the bootstrap bands `VARirband` produces.
- `LPmodel` projects only `ENDO`'s **first column**.
- `zscore` (used for `impact = 0`) and `norminv` are Statistics Toolbox
  functions — two more toolbox dependencies beyond those catalogued in the
  brief. Mapped to `(s .- mean(s)) ./ std(s)` and `quantile(Normal(), ·)`.

**Not checked:** `LPmodel` with multi-column `S`.
