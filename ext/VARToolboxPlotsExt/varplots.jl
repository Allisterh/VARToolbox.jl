# =======================================================================
# varplots.jl — the six VAR/SR plot functions
#
# Sources: VAR/VARirplot.m, VAR/SRirplot.m, VAR/VARvdplot.m, VAR/SRvdplot.m,
# VAR/VARhdplot.m, VAR/SRhdplot.m (v3dot1), with the `quality==2` save
# branches restored from their v3dot0 counterparts (v3dot1 dropped them
# while defaulting `quality = 2`, so upstream saves nothing by default).
# All other deviations are marked `DEVIATION` inline and listed in
# PORTING_LOG.md.
# =======================================================================

# MATLAB `round(sqrt(n))` rounds ties away from zero.
_grid_dims(nvars) = (_round_half_away(sqrt(nvars)), ceil(Int, sqrt(nvars)))

# `char('a'+(1:n)-1)` then indexing — the suptitle letter for panel jj.
_alphabet(jj) = string(Char('a' + jj - 1))

# A (row×col) panel grid, sized in cm, unused trailing panels blanked.
function _panel_grid(row, col, npanels, figsize)
    p = plot(layout = (row, col), size = FigSize(figsize[1], figsize[2]))
    for ii in 1:(row * col)
        ii <= npanels ? _matlab_axes!(p[ii]) : plot!(p[ii]; framestyle = :none)
    end
    return p
end

# Compose a figure with a bottom legend strip (LegSubplot equivalent —
# MATLAB shrinks the axes and floats one centred box-less legend; Plots
# needs a dedicated subplot row).
function _with_legend_row(pmain, entries; size = pmain.attr[:size])
    pleg = plot(framestyle = :none)
    _legend_row!(pleg, entries)
    return plot(pmain, pleg; layout = Plots.grid(2, 1, heights = [0.9, 0.1]),
                size = size)
end

# In-axes legend for the `quality == 0` branches (MATLAB `legend(...)` on the
# current — i.e. last — axes, default 'northeast' placement).
function _inline_legend!(sp, entries)
    plot!(sp; legend = :topright, legendfontsize = 8,
          legend_background_color = :white)
    for (lab, spec) in entries
        if spec isa RGB
            scatter!(sp, [NaN]; marker = :square, markersize = 6,
                     markercolor = spec, markerstrokecolor = spec, label = lab)
        else
            plot!(sp, [NaN]; linecolor = spec.color, linewidth = spec.lw,
                  label = lab)
        end
    end
    return sp
end

_area_legend_entries(snames) = [(s, cmap(k)) for (k, s) in enumerate(snames)]

# ---- VARirplot / SRirplot ---------------------------------------------

# Shared body: the two differ only in name checks, the number of shock
# figures, the colour/marker of the (band-hidden) first median line, and the
# suptitle label source.
function _irplot(IR, VARopt, INF, SUP; sr::Bool)
    vnames = VARopt.vnames
    isempty(vnames) &&
        error("You need to add label for endogenous variables in VARopt")
    if sr
        isempty(VARopt.snames) &&
            error("You need to add label for shocks in VARopt")
        snames = VARopt.snames
        length(snames) > length(vnames) &&
            error("You have more shock labels than variables")
    else
        snames = isempty(VARopt.snames) ? vnames : VARopt.snames
    end

    filename = VARopt.figname * "IR_"
    quality = VARopt.quality
    suptitle = VARopt.suptitle
    pick = VARopt.pick

    nsteps, nvars, _ = size(IR)
    # SRirplot.m:52 counts shocks off snames, not size(IR,3).
    nshocks = sr ? length(snames) : size(IR, 3)

    if pick < 0 || pick > nvars
        error("The selected shock is non valid")
    elseif pick == 0
        pick = 1
    else
        nshocks = pick
    end

    row, col = _grid_dims(nvars)
    steps = 1:nsteps

    plots = Plots.Plot[]
    for jj in pick:nshocks
        # DEVIATION: fresh figure per shock. v3dot1 comments out `clf('reset')`
        # in VARirplot so successive shocks overplot into one axes; v3dot0
        # (and the toolbox's committed reference figures) clear between
        # shocks. SRirplot's clf is active even in v3dot1.
        p = _panel_grid(row, col, nvars, VARopt.FigSize)
        for ii in 1:nvars
            sp = p[ii]
            # First median line, drawn under the swathe exactly as upstream
            # (black '*' in VARirplot.m:82, cmap(1) 'x' in SRirplot.m:83);
            # PlotSwathe re-draws its own navy line on top, hiding this one.
            firstcol = sr ? cmap(1) : RGB(0.0, 0.0, 0.0)
            plot!(sp, steps, IR[:, ii, jj]; linecolor = firstcol,
                  linewidth = 2, marker = sr ? :x : :star6, markersize = 4,
                  markercolor = firstcol, markerstrokecolor = firstcol,
                  label = "")
            if INF !== nothing && SUP !== nothing
                PlotSwathe!(sp, IR[:, ii, jj],
                            hcat(INF[:, ii, jj], SUP[:, ii, jj]);
                            marker = :star6)
            end
            plot!(sp, steps, zeros(nsteps); linecolor = :black,
                  linestyle = :dash, linewidth = 0.5, label = "")
            plot!(sp; xlims = (1, nsteps),
                  title = vnames[ii] * " to " * snames[jj])
        end
        FigName = filename * string(jj)
        if quality == 1 && suptitle == 1
            # VARirplot.m:96 labels the shock with vnames{jj} — an upstream
            # slip replicated verbatim; SRirplot.m:97 uses snames{jj}.
            lab = sr ? snames[jj] : vnames[jj]
            plot!(p; plot_title = _alphabet(jj) * ") IR to a shock to " * lab,
                  plot_titlefontsize = 14)
        end
        SaveFigure(p, FigName, quality)
        push!(plots, p)
    end
    return plots
end

VARirplot(IR::AbstractArray{<:Real,3}, VARopt::VAROptions,
          INF = nothing, SUP = nothing) = _irplot(IR, VARopt, INF, SUP; sr = false)

SRirplot(IR::AbstractArray{<:Real,3}, VARopt::VAROptions,
         INF = nothing, SUP = nothing) = _irplot(IR, VARopt, INF, SUP; sr = true)

# ---- VARvdplot / SRvdplot ---------------------------------------------

function _vdplot(VD, VARopt; sr::Bool)
    vnames = VARopt.vnames
    isempty(vnames) &&
        error("You need to add label for endogenous variables in VARopt")
    if sr
        isempty(VARopt.snames) &&
            error("You need to add label for shocks in VARopt")
        snames = VARopt.snames
    else
        snames = isempty(VARopt.snames) ? vnames : VARopt.snames
    end

    filename = VARopt.figname * "VD"
    quality = VARopt.quality
    suptitle = VARopt.suptitle
    pick = VARopt.pick

    # Upstream destructures VD's (nsteps × nshocks × nvar) layout under the
    # names [nsteps, nvars, nshocks] (VARvdplot.m:52) — replicated verbatim,
    # including the resulting role swap; harmless for the square VD every
    # VARvd call produces.
    nsteps, nvars, nshocks = size(VD)

    if pick < 0 || pick > nvars
        error("The selected variable is non valid")
    elseif pick == 0
        pick = 1
    else
        nvars = pick
    end

    row, col = _grid_dims(nshocks)

    p = _panel_grid(row, col, nvars, VARopt.FigSize)
    for ii in 1:nvars
        sp = p[ii]
        AreaPlot!(sp, VD[:, :, ii])
        plot!(sp; xlims = (1, nsteps), ylims = (0, 100), ywiden = false,
              title = vnames[ii])
    end

    entries = _area_legend_entries(snames)
    FigName = filename
    if quality == 1 || quality == 2
        if suptitle == 1
            # VARvdplot.m:90 reads the loop-terminal `ii`; SRvdplot.m prints a
            # fixed title.
            st = sr ? "Variance Decomposition" :
                      _alphabet(nvars) * ") VD of " * vnames[nvars]
            plot!(p; plot_title = st, plot_titlefontsize = 14)
        end
        p = _with_legend_row(p, entries;
                             size = FigSize(VARopt.FigSize[1], VARopt.FigSize[2]))
    else
        _inline_legend!(p[nvars], entries)
    end
    SaveFigure(p, FigName, quality)
    return p
end

VARvdplot(VD::AbstractArray{<:Real,3}, VARopt::VAROptions) =
    _vdplot(VD, VARopt; sr = false)

SRvdplot(VD::AbstractArray{<:Real,3}, VARopt::VAROptions) =
    _vdplot(VD, VARopt; sr = true)

# ---- VARhdplot / SRhdplot ---------------------------------------------

function _hdplot(HD::HDOutput, VARopt; sr::Bool)
    vnames = VARopt.vnames
    isempty(vnames) &&
        error("You need to add label for endogenous variables in VARopt")
    if sr
        isempty(VARopt.snames) &&
            error("You need to add label for shocks in VARopt")
        snames = VARopt.snames
    else
        snames = isempty(VARopt.snames) ? vnames : VARopt.snames
    end

    filename = VARopt.figname * "HD_"
    quality = VARopt.quality
    suptitle = VARopt.suptitle
    pick = VARopt.pick

    # Same verbatim role-swapped destructure as the VD plotters
    # (VARhdplot.m:46); nsteps here is nobs + nlag.
    nsteps, nvars, nshocks = size(HD.shock)

    if pick < 0 || pick > nvars
        error("The selected shock is non valid")
    elseif pick == 0
        pick = 1
    else
        nshocks = pick
    end

    plots = Plots.Plot[]
    for ii in pick:nvars
        p = plot(size = FigSize(VARopt.FigSize[1], VARopt.FigSize[2]))
        _matlab_axes!(p)
        AreaPlot!(p, HD.shock[:, :, ii])
        plot!(p, 1:nsteps, vec(sum(HD.shock[:, :, ii], dims = 2));
              linecolor = :black, linewidth = 2, label = "")
        VARopt.firstdate !== nothing &&
            DatesPlot!(p, VARopt.firstdate, nsteps, 8, VARopt.frequency)
        plot!(p; xlims = (1, nsteps), title = vnames[ii])

        entries = vcat(_area_legend_entries(snames),
                       [("Data", (color = RGB(0.0, 0.0, 0.0), lw = 2))])
        FigName = filename * string(ii)
        if quality == 1 || quality == 2
            if suptitle == 1
                # DEVIATION: upstream reads an undefined loop index `jj` here
                # (VARhdplot.m:75) and always crashes; the evident intent, the
                # current variable, is substituted.
                plot!(p; plot_title = _alphabet(ii) * ") HD of " * vnames[ii],
                      plot_titlefontsize = 14)
            end
            p = _with_legend_row(p, entries;
                                 size = FigSize(VARopt.FigSize[1], VARopt.FigSize[2]))
        else
            # DEVIATION: VARhdplot.m:81 references an undefined `vsnames` and
            # always crashes; `snames` is the evident intent. SRhdplot.m:83
            # genuinely (and buggily) labels with `vnames` — replicated.
            leg = sr ? vcat(_area_legend_entries(vnames),
                            [("Data", (color = RGB(0.0, 0.0, 0.0), lw = 2))]) :
                       entries
            _inline_legend!(p, leg)
        end
        SaveFigure(p, FigName, quality)
        push!(plots, p)
    end
    return plots
end

VARhdplot(HD::HDOutput, VARopt::VAROptions) = _hdplot(HD, VARopt; sr = false)

SRhdplot(HD::HDOutput, VARopt::VAROptions) = _hdplot(HD, VARopt; sr = true)
