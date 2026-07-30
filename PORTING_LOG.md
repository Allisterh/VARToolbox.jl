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
