# =======================================================================
# primitives.jl — ports of the Figure/ layer building blocks
#
# Sources: Figure/cmap.m, Figure/PlotSwathe.m (+ PlotSwatheOption.m),
# Figure/AreaPlot.m, Figure/DatesPlot.m, Figure/FigSize.m,
# Figure/SaveFigure.m, Figure/LegOption.m + LegSubplot.m.
# =======================================================================

# ---- cmap -------------------------------------------------------------

# The nine fixed rows of Figure/cmap.m, exact.
const _CMAP = [
    RGB(0.000, 0.447, 0.741),        # 1 blue
    RGB(0.929, 0.694, 0.125),        # 2 yellow
    RGB(0.494, 0.184, 0.556),        # 3 purple
    RGB(0.301, 0.745, 0.933),        # 4 light blue
    RGB(0.850, 0.325, 0.098),        # 5 red
    RGB(0.466, 0.674, 0.188),        # 6 green
    RGB(0.635, 0.078, 0.184),        # 7 dark red
    RGB(0 / 255, 48 / 255, 143 / 255),   # 8 air force blue
    RGB(209 / 255, 96 / 255, 1 / 255),   # 9 orange
]

function cmap(n::Integer)
    1 <= n <= 9 && return _CMAP[n]
    error("cmap($n): only the 9 fixed colours are ported; upstream falls " *
          "back to MATLAB's parula colormap for n > 9 (see PORTING_LOG.md). " *
          "Supply your own colour instead.")
end

# MATLAB PlotSwatheOption.m defaults.
const _SWATHE_LINECOL = RGB(13 / 255, 54 / 255, 84 / 255)      # dark navy
const _SWATHE_COL = RGB(138 / 255, 178 / 255, 212 / 255)       # steel blue

# PlotSwathe.m:98 — band edge lines are the swathe colour pulled 55% toward
# white.
_lighten(c::RGB, fac) = RGB(c.r + (1 - c.r) * fac,
                            c.g + (1 - c.g) * fac,
                            c.b + (1 - c.b) * fac)

# ---- MATLAB-look axes defaults ---------------------------------------

# MATLAB axes: boxed frame, no grid, Helvetica 10pt, inward ticks, bold
# 10pt titles (every toolbox call passes 'FontWeight','bold','FontSize',10).
# GR's default font is already Helvetica-like; forcing fontfamily="Helvetica"
# routes GR through a TTF whose Unicode minus (U+2212) is missing, turning
# every negative tick label into "?" — so only the title family is pinned
# (titles carry no minus signs).
function _matlab_axes!(sp)
    plot!(sp; framestyle = :box, grid = false, tick_direction = :in,
          tickfontsize = 10,
          titlefontsize = 10, titlefontfamily = "Helvetica Bold",
          legend = false, xwiden = false)
    return sp
end

# ---- FigSize ----------------------------------------------------------

# GR maps `size` pixels 1:1 to PDF points (verified: MediaBox == size), so a
# page physically matching MATLAB's cm window needs px = cm * 72 / 2.54.
FigSize(x_cm::Real = 26, y_cm::Real = 24) =
    (round(Int, x_cm * 72 / 2.54), round(Int, y_cm * 72 / 2.54))

# ---- SaveFigure -------------------------------------------------------

# Upstream quality 0/1/2 dispatch collapses to one native vector export; the
# argument is kept for call-site compatibility. See PORTING_LOG.md.
function SaveFigure(p::Plots.Plot, figname::AbstractString, quality::Integer = 2)
    savefig(p, figname * ".pdf")
    return nothing
end

# ---- PlotSwathe -------------------------------------------------------

function PlotSwathe!(sp, bar::AbstractVector, swathe;
                     linecol = _SWATHE_LINECOL,
                     swathecol = _SWATHE_COL,
                     linewidth::Real = 2,
                     marker = :none,
                     transp::Bool = false,
                     alpha::Real = 1.0,
                     swatheonly::Bool = false,
                     label::AbstractString = "")
    nobs = length(bar)
    barv = vec(collect(float.(bar)))

    # PlotSwathe.m:71-91 — vector swathe is symmetric around `bar`; a matrix
    # is a pair of band limits. MATLAB flips an n×2 [INF SUP] to 2×n and then
    # names row 1 (INF) "upp" and row 2 (SUP) "low"; the swap is harmless for
    # a filled band and is replicated by simply filling between the two.
    if swathe isa AbstractVector
        sw = vec(collect(float.(swathe)))
        low = barv .- sw
        upp = barv .+ sw
    else
        swm = collect(float.(swathe))
        size(swm, 1) > 2 && (swm = permutedims(swm))
        upp = vec(swm[1, :])
        low = vec(swm[2, :])
    end

    xaxis = 1:nobs

    # PlotSwathe.m:98-110 — edge colour is lightened; with the default
    # transp=0, alpha=1 the patch colour is exactly swathecol, fully opaque.
    edgecolor = _lighten(swathecol, 0.55)
    if transp
        facealpha = alpha
        patchcolor = swathecol
    else
        facealpha = 1.0
        patchcolor = _lighten(swathecol, 1 - alpha)
    end

    # MATLAB strips NaN vertices from the patch polygon (PlotSwathe.m:126-127),
    # joining across gaps; masking both limits reproduces that.
    mask = .!(isnan.(low) .| isnan.(upp))

    # Net MATLAB z-order (line drawn, deleted, re-drawn on top): patch, then
    # edge lines, then the bar line.
    plot!(sp, xaxis[mask], low[mask]; fillrange = upp[mask],
          fillcolor = patchcolor, fillalpha = facealpha,
          linealpha = 0, label = "")
    plot!(sp, xaxis, low; linecolor = edgecolor, linewidth = 0.5, label = "")
    plot!(sp, xaxis, upp; linecolor = edgecolor, linewidth = 0.5, label = "")
    if !swatheonly
        plot!(sp, xaxis, barv; linecolor = linecol, linewidth = linewidth,
              marker = marker, markersize = 4, markercolor = linecol,
              markerstrokecolor = linecol, label = label)
    end
    return sp
end

PlotSwathe!(bar::AbstractVector, swathe; kwargs...) =
    PlotSwathe!(Plots.current(), bar, swathe; kwargs...)

function PlotSwathe(bar::AbstractVector, swathe; kwargs...)
    p = plot()
    _matlab_axes!(p)
    PlotSwathe!(p, bar, swathe; kwargs...)
    return p
end

# ---- AreaPlot ---------------------------------------------------------

function AreaPlot!(sp, X::AbstractMatrix; labels = nothing)
    Xf = collect(float.(X))
    n, nser = size(Xf)
    x = 1:n
    # AreaPlot.m:26-36 — positive part stacked above zero, negative below,
    # both in cmap order with no edge line. MATLAB's area() stacks columns
    # cumulatively; NaN rows (e.g. VARhd's first nlag observations) stay NaN
    # through the row-wise cumsum and render as gaps, as in MATLAB.
    for (part, labelled) in ((Xf .* (Xf .> 0), true), (Xf .* (Xf .< 0), false))
        C = cumsum(part, dims = 2)
        for k in 1:nser
            lower = k == 1 ? zeros(n) : C[:, k - 1]
            lab = (labelled && labels !== nothing) ? labels[k] : ""
            plot!(sp, x, C[:, k]; fillrange = lower, fillcolor = cmap(k),
                  fillalpha = 1, linealpha = 0, label = lab)
        end
    end
    return sp
end

# ---- DatesPlot --------------------------------------------------------

# MATLAB round() rounds ties away from zero; Julia's default rounds to even.
_round_half_away(x) = round(Int, x, RoundNearestTiesAway)

"""
    _date_ticks(fo, nobs, nticks, frequency) -> (positions, labels, xlims)

Pure transliteration of Figure/DatesPlot.m: tick positions on a `1:nobs`
axis, year labels, and the `(0, nobs+1)` x-limits it imposes. Kept separate
from the mutating [`DatesPlot!`](@ref) so it can be unit-tested against
MATLAB output.

The monthly branch follows VAR Toolbox 4.0, which replaced 3.1's literal
`0.0833` step with exact `1/12` arithmetic and `round`-based month parsing
(`ChangeLog.md`, v4.0 release, `Figure/DatesPlot.m`). Under 3.1 the
accumulated `0.0833` drift made every monthly tick position fractional
(12.0048, 24.0048, ...) and MATLAB floored them, with a console warning,
when indexing the label vector.
"""
function _date_ticks(fo::Real, nobs::Integer, nticks::Real, frequency::AbstractString = "q")
    frequency in ("y", "q", "m") ||
        error("Frequency should be either 'y', 'q', or 'm'")

    if frequency == "y"
        lo = fo + nobs - 1
        labyear = collect(fo:lo)
        nobslabyear = length(labyear)
        tick = _round_half_away(nobs / nticks)
        tick_max = floor(Int, nticks) * tick
        lo_new = lo
        for _ in 1:(tick_max - nobslabyear)
            lo_new += 1
        end
        labyear = collect(fo:lo_new)
        idx = tick:tick:tick_max
        positions = collect(Float64, idx)
        labels = labyear[collect(idx)]
    elseif frequency == "q"
        step = 0.25
        per = 4
        # DatesPlot.m accumulates the last observation with a loop; replicate
        # so float error matches.
        lo = fo
        for _ in 1:(nobs - 1)
            lo += step
        end
        fo_year = floor(fo)
        fo_diff = abs(fo - fo_year - 0.75) / step + 1
        labyear = floor.(collect(fo:step:lo))
        nobslabyear = length(labyear)
        lo_year = (lo - floor(lo) == 0.75) ? floor(lo) + 1 : floor(lo)
        nobs_year = lo_year - fo_year + 1
        tick = _round_half_away(nobs_year / nticks)
        if floor(nticks) == nticks
            tick_max = floor(Int, nticks) * (tick * per)
        else
            tick_max = floor(Int, nticks) * (tick * per) + tick * per
        end
        lo_new = lo
        for _ in 1:(tick_max - nobslabyear)
            lo_new += step
        end
        labyear_new = floor.(collect((lo + step):step:lo_new))
        labyear = vcat(labyear, labyear_new)
        positions = collect(fo_diff:(tick * per):tick_max)
        labels = labyear[_round_half_away.(positions)]
    else
        # Monthly, v4.0: exact 1/12 arithmetic throughout.
        per = 12
        lo = fo + (nobs - 1) / 12
        fo_year = floor(fo)
        fo_diff = 12 - _round_half_away(12 * (fo - fo_year))
        lo_year = _round_half_away(12 * (lo - floor(lo))) == 11 ? floor(lo) + 1 : floor(lo)
        labyear = floor.(collect(fo:(1/12):lo))
        nobslabyear = length(labyear)
        nobs_year = lo_year - fo_year + 1
        tick = _round_half_away(nobs_year / nticks)
        if floor(nticks) == nticks
            tick_max = floor(Int, nticks) * (tick * per)
        else
            tick_max = floor(Int, nticks) * (tick * per) + tick * per
        end
        lo_new = lo + (tick_max - nobslabyear) / 12
        labyear_new = floor.(collect((lo + 1/12):(1/12):lo_new))
        labyear = vcat(labyear, labyear_new)
        positions = collect(Float64, fo_diff:(tick * per):tick_max)
        labels = labyear[Int.(positions)]
    end

    return positions, string.(Int.(labels)), (0.0, nobs + 1.0)
end

function DatesPlot!(sp, fo::Real, nobs::Integer, nticks::Real,
                    frequency::AbstractString = "q")
    positions, labels, lims = _date_ticks(fo, nobs, nticks, frequency)
    plot!(sp; xticks = (positions, labels), xlims = lims)
    return sp
end

# ---- Bottom legend strip (LegOption/LegSubplot equivalent) ------------

# LegSubplot.m places one centred, box-less horizontal legend along the
# bottom of the figure. Plots has no figure-level legend, so figures that
# need one are built with a dedicated frameless subplot row; this fills it
# with invisible proxy series carrying the right colours and labels.
# `entries` is a vector of (label, spec) where spec is either an RGB (area
# swatch, drawn as a filled square marker) or a NamedTuple (line entry).
function _legend_row!(sp, entries; columns = length(entries))
    plot!(sp; framestyle = :none, legend = :top, legend_columns = columns,
          legend_background_color = :transparent,
          legend_foreground_color = :transparent,
          legendfontsize = 9)
    for (lab, spec) in entries
        if spec isa RGB
            scatter!(sp, [NaN]; marker = :square, markersize = 7,
                     markercolor = spec, markerstrokecolor = spec, label = lab)
        else
            plot!(sp, [NaN]; linecolor = spec.color, linewidth = spec.lw,
                  label = lab)
        end
    end
    return sp
end
