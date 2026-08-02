# =======================================================================
# plotstubs.jl — function stubs for the plotting layer
#
# The plotting layer is a package extension (ext/VARToolboxPlotsExt) with a
# weak dependency on Plots.jl, so core estimation never drags in a plotting
# stack. These stubs own the names and docstrings; the methods are added by
# the extension, which loads automatically once the user does `using Plots`.
#
# Ported from VAR Toolbox 3.1 (Ambrogio Cesa-Bianchi), files VAR/VARirplot.m,
# VAR/VARvdplot.m, VAR/VARhdplot.m, VAR/SRirplot.m, VAR/SRvdplot.m,
# VAR/SRhdplot.m and Figure/{PlotSwathe,AreaPlot,cmap,DatesPlot,FigSize,
# SaveFigure}.m. GPL-3.0, see PROVENANCE.md.
# =======================================================================

# Internal: throw a uniform, actionable error from zero-method stubs when a
# user calls a plotting function without having loaded Plots.
function _plots_not_loaded(name::Symbol)
    error("$name is part of the plotting layer, which activates when Plots " *
          "is loaded. Run `using Plots` first (extension VARToolboxPlotsExt).")
end

"""
    VARirplot(IR, VARopt, INF=nothing, SUP=nothing) -> Vector{Plots.Plot}

Plot impulse responses computed with [`VARir`](@ref) — one figure **per
shock**, each saved as `[VARopt.figname]IR_<shock#>.pdf`. Panels are laid out
`round(sqrt(nvar)) × ceil(sqrt(nvar))`, one per response variable, each with
the median line (dark navy, `*`-style marker), the `INF`/`SUP` swathe when
given, and a dashed zero line. Titles are `"<vname> to <sname>"`.

Reads off `VARopt`: `vnames` (required), `snames` (falls back to `vnames`),
`figname`, `FigSize` (cm), `quality`, `suptitle`, `pick`.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note "`pick` selects a single shock here"
    `pick = k > 0` plots **only** shock `k`. This is inconsistent with
    [`VARvdplot`](@ref) (variables `1:k`) and [`VARhdplot`](@ref) (variables
    `k:nvar`); all three upstream meanings are replicated verbatim.

!!! note "Deviations from v3dot1 (documented in PORTING_LOG.md)"
    v3dot1 dropped the `quality==2` save branch (the `VARoption` default), so
    upstream saves nothing by default; this port restores v3dot0's behaviour
    of writing a vector PDF. v3dot1 also comments out `clf('reset')`, so
    successive shocks overplot into one axes upstream; this port gives each
    shock a fresh figure, which is what v3dot0 did and what the toolbox's own
    committed reference figures show.
"""
function VARirplot end

"""
    VARvdplot(VD, VARopt) -> Plots.Plot

Plot the forecast error variance decomposition from [`VARvd`](@ref) as
stacked areas (one panel per variable, shocks in [`cmap`](@ref) colours,
y-axis 0–100), saved as `[VARopt.figname]VD.pdf` with a horizontal legend
strip of `snames` along the bottom.

Reads off `VARopt`: `vnames` (required), `snames` (falls back to `vnames`),
`figname`, `FigSize`, `quality`, `suptitle`, `pick` (`pick = k > 0` plots
variables `1:k`).

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note
    Upstream destructures `size(VD)` as `[nsteps, nvars, nshocks]` even
    though `VD` is laid out `(nsteps × nshocks × nvar)` — harmless for the
    square `VD` every `VARvd` call produces, and replicated verbatim. The
    missing v3dot1 `quality==2` save branch is restored from v3dot0.
"""
function VARvdplot end

"""
    VARhdplot(HD, VARopt) -> Vector{Plots.Plot}

Plot the historical decomposition from [`VARhd`](@ref) — one full-size figure
per variable, saved as `[VARopt.figname]HD_<var#>.pdf`: stacked shock
contributions ([`AreaPlot!`](@ref)) with the actual (demeaned) data overlaid
as a black line (the sum of contributions), and year tick labels via
[`DatesPlot!`](@ref) when `VARopt.firstdate` is set.

Reads off `VARopt`: `vnames` (required), `snames` (falls back to `vnames`),
`figname`, `FigSize`, `quality`, `suptitle`, `pick` (`pick = k > 0` plots
variables `k:nvar`), `firstdate`, `frequency`.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note "Two upstream always-crash branches are repaired"
    Upstream's `quality==0` legend references an undefined `vsnames` and its
    `suptitle==1` branch an undefined loop index; both crash MATLAB
    unconditionally when reached, so this port substitutes the evident
    intent (`snames`, current variable) — the only kind of deviation the
    port permits itself. See PORTING_LOG.md.
"""
function VARhdplot end

"""
    SRirplot(IR, VARopt, INF=nothing, SUP=nothing) -> Vector{Plots.Plot}

[`VARirplot`](@ref) for sign-restriction output from [`SR`](@ref): pass
`SRout.IR` (or `SRout.IRmed`) with `SRout.IRinf`/`SRout.IRsup`. Differences
from `VARirplot`: both `vnames` **and** `snames` are required, and the number
of shock figures is `length(snames)` — not `size(IR, 3)` — so a
single-restricted-shock run yields a single figure.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function SRirplot end

"""
    SRvdplot(VD, VARopt) -> Plots.Plot

[`VARvdplot`](@ref) for sign-restriction output from [`SR`](@ref); both
`vnames` and `snames` are required. Saved as `[VARopt.figname]VD.pdf`.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function SRvdplot end

"""
    SRhdplot(HD, VARopt) -> Vector{Plots.Plot}

[`VARhdplot`](@ref) for sign-restriction output from [`SR`](@ref); both
`vnames` and `snames` are required.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function SRhdplot end

"""
    PlotSwathe(bar, swathe; kwargs...)  -> Plots.Plot
    PlotSwathe!([sp], bar, swathe; kwargs...) -> sp

Line with a shaded swathe (`Figure/PlotSwathe.m`). `swathe` is either a
vector (symmetric band `bar ± swathe`) or an `n × 2` matrix `[INF SUP]` of
band limits. Draw order matches MATLAB: opaque band fill, lightened band
edge lines, then the `bar` line on top.

Keywords (MATLAB `PlotSwatheOption` defaults): `linecol` (dark navy
`RGB(13/255, 54/255, 84/255)`), `swathecol` (steel blue
`RGB(138/255, 178/255, 212/255)`), `linewidth = 2`, `marker = :none`,
`transp = false`, `alpha = 1.0`, `swatheonly = false`, `label = ""`.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note
    The band is **opaque** by default, exactly as upstream: MATLAB's
    `VARirplot` sets the misspelled field `SwatheOpt.trans` (the real field
    is `transp`), so transparency is silently never enabled there.
"""
function PlotSwathe end

@doc (@doc PlotSwathe)
function PlotSwathe! end

"""
    AreaPlot!(sp, X; labels=nothing) -> sp

Stacked area chart with positive values stacked above zero and negative
values stacked below (`Figure/AreaPlot.m`), series `k` filled in `cmap(k)`
with no edge line. `labels` (a vector of strings) labels the positive-stack
series for the legend, matching MATLAB's use of the first handle row.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function AreaPlot! end

"""
    cmap(n) -> RGB

The toolbox's ordinal colour palette (`Figure/cmap.m`): nine fixed colours
(1 blue, 2 yellow, 3 purple, 4 light blue, 5 red, 6 green, 7 dark red,
8 air-force blue, 9 orange).

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note
    Upstream falls back to rows of MATLAB's `parula` colormap for `n > 9`;
    that tail is not ported (no toolbox figure reaches it) and `n > 9`
    throws. See PORTING_LOG.md.
"""
function cmap end

"""
    DatesPlot!(sp, fo, nobs, nticks, frequency="q") -> sp

Year tick labels on the x-axis of a chart indexed `1:nobs`
(`Figure/DatesPlot.m`). `fo` is the numeric first observation under the
toolbox convention `1987.00 = 1987Q1` (quarterly step `0.25`, monthly step
`0.0833`, annual integers). Sets `xlims = (0, nobs+1)` exactly as upstream.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function DatesPlot! end

"""
    FigSize(x_cm=26, y_cm=24) -> (w_px, h_px)

MATLAB's `Figure/FigSize.m` sizes figure windows in centimetres; Plots sizes
figures in pixels that GR maps 1:1 to PDF points. This helper converts:
`FigSize(26, 12)` returns the pixel tuple whose saved PDF page measures
26 cm × 12 cm, for use as `plot(...; size = FigSize(26, 12))`.

Requires `using Plots` (package extension `VARToolboxPlotsExt`).
"""
function FigSize end

"""
    SaveFigure(p, figname, quality=2)

Write plot `p` to `figname.pdf` (`Figure/SaveFigure.m`).

Requires `using Plots` (package extension `VARToolboxPlotsExt`).

!!! note "All quality levels write the same vector PDF"
    Upstream dispatches on `quality`: `0` → rasterised `print -r100`, `1` →
    Ghostscript `export_fig` (not shipped in v3dot1 — always crashes), `2` →
    `exportgraphics`. Julia's native vector export supersedes all three, so
    the argument is accepted for signature compatibility and ignored.
    See PORTING_LOG.md.
"""
function SaveFigure end
