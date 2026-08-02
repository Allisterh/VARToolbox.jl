# =======================================================================
# figures.jl
#
# Reproduce, from the Julia port, every figure the MATLAB toolbox produces
# for the four published replications (v3dot1/Replic/*/GO_*.m), plus the
# Primer §9 historical-decomposition figures — the only place upstream
# exercises VARhdplot and date axes.
#
# The bootstrap/sign-restriction draws are replayed from the committed
# golden-master fixture test/fixtures/phase6.mat (recorded from MATLAB
# R2026a), so every plotted curve is numerically identical to a MATLAB run
# with those draws — any visual difference against the MATLAB-rendered
# reference is a rendering difference, not a numerical one.
#
# Requires the plotting extension: `using Plots` must precede inclusion
# (see examples/Project.toml). Regenerate everything with:
#
#     julia --project=julia/VARToolbox/examples -e '
#         include("julia/VARToolbox/examples/figures.jl")
#         generate_all_figures()'
#
# One deliberate deviation, documented in figures/README.md: the Uhlig
# figure uses the fixture's 100 accepted draws, not GO_Uhlig2005.m's 500 —
# the committed fixture was capped for size. Bands are slightly less smooth;
# the identification content is unchanged.
# =======================================================================

using VARToolbox
using Plots
using MAT
using Random

isdefined(@__MODULE__, :load_replication_data) ||
    include(joinpath(@__DIR__, "replications.jl"))

const FIGURES_DIR = joinpath(@__DIR__, "figures")
const FIXTURES_DIR = normpath(joinpath(@__DIR__, "..", "test", "fixtures"))
const PRIMER_DIR = normpath(joinpath(@__DIR__, "..", "data", "Primer"))

_pool(x) = [Matrix{Float64}(a) for a in vec(x)]
_postpool(x) = [(Matrix{Float64}(p[1]), Matrix{Float64}(p[2])) for p in vec(x)]

# MATLAB `axis tight`: both axes shrink to the exact data extent. The zero
# line is itself plotted data, so the y-range always includes 0. Ticks are
# thinned explicitly: MATLAB keeps labels short with an axis multiplier
# (x10^-3) that GR does not support, so unthinned decimal labels overlap in
# short panels.
function _axis_tight!(sp, nsteps, series...)
    lo = min(0.0, minimum(minimum.(series)))
    hi = max(0.0, maximum(maximum.(series)))
    ticks = Plots.optimize_ticks(lo, hi; k_min = 3, k_max = 5)[1]
    plot!(sp; xlims = (1, nsteps), ylims = (lo, hi),
          xwiden = false, ywiden = false, yticks = ticks)
    return sp
end

_ext() = Base.get_extension(VARToolbox, :VARToolboxPlotsExt)

"""
    primer_sw_section9(outdir) -> Vector{Plots.Plot}

Primer §9.6: historical decomposition of the Stock-Watson VAR — the Primer
orders the variables (inflation, unemployment, fed funds), unlike
`GO_SW2001.m`'s spreadsheet order. Produces `SW_HD_1..3.pdf`, the only
upstream figures that exercise `VARhdplot` and `DatesPlot`. Deterministic —
no draws involved.
"""
function primer_sw_section9(outdir::AbstractString)
    path = joinpath(PRIMER_DIR, "SW2001_Data.xlsx")
    sheet = XLSX.readxlsx(path)["Sheet1"]
    raw = sheet[:]
    vnames = String[string(x) for x in raw[2, 2:end]]
    dates = String[string(x) for x in raw[3:end, 1]]
    body = raw[3:end, 2:end]
    data = Matrix{Float64}(undef, size(body)...)
    for i in eachindex(body)
        v = body[i]
        data[i] = (v === missing) ? NaN : Float64(v)
    end
    data = Num2NaN(data)

    # Date2Num: '1960q1' -> 1960 + 1/4 - 0.25 = 1960.00.
    m = match(r"^(\d{4})[Qq](\d)$", dates[1])
    firstdate = parse(Float64, m.captures[1]) + (parse(Int, m.captures[2]) - 1) * 0.25

    X = hcat((col(data, vnames, n) for n in ("infl", "unemp", "ff"))...)

    VAR, VARopt = VARmodel(X, 4, 1)
    VARopt = @set VARopt.vnames = ["Inflation", "Unemployment", "Fed Funds"]
    VARopt = @set VARopt.nsteps = 24
    VARopt = @set VARopt.ident = "short"
    # The Primer's snames are TeX ('\epsilon^{1}', ...); GR renders plain
    # strings literally, so the closest Unicode is used instead.
    VARopt = @set VARopt.snames = ["ε¹", "ε²", "ε^MonPol"]
    VARopt = @set VARopt.FigSize = [26, 12]
    VARopt = @set VARopt.firstdate = firstdate
    VARopt = @set VARopt.frequency = "q"
    VARopt = @set VARopt.figname = joinpath(outdir, "SW_")

    IR, VAR = VARir(VAR, VARopt)
    HD, VAR = VARhd(VAR, VARopt)
    return VARhdplot(HD, VARopt)
end

"""
    generate_all_figures(outdir=FIGURES_DIR; fixtures=FIXTURES_DIR)
        -> Dict{String,Any}

Generate all 12 replication figures as PDFs in `outdir`, returning the
`Plots.Plot` objects keyed by output stem. See the file header for what
"replication" means numerically.
"""
function generate_all_figures(outdir::AbstractString = FIGURES_DIR;
                              fixtures::AbstractString = FIXTURES_DIR)
    mkpath(outdir)
    R = matread(joinpath(fixtures, "phase6.mat"))
    out = Dict{String,Any}()
    E = _ext()
    E === nothing && error("The plotting extension is not loaded — `using Plots` first.")

    # ---- Stock & Watson (2001): VARirplot + VARvdplot -----------------
    # GO_SW2001.m: vnames = long names from the sheet, FigSize [26,12],
    # snames left empty (falls back to vnames), plots IRbar with the bands.
    _, _, sw_long, _ = load_replication_data(joinpath("SW2001", "SW2001_Data.xlsx"))
    sw = stock_watson_2001(ndraws = Int(R["sw_ndraws"]),
                           band_src = DrawSource(uniforms = _pool(R["sw_irband_draws"])),
                           vd_src = DrawSource(uniforms = _pool(R["sw_vdband_draws"])))
    opt = sw.VARopt
    opt = @set opt.vnames = sw_long
    opt = @set opt.FigSize = [26, 12]
    opt = @set opt.figname = joinpath(outdir, "SW2001_")
    ps = VARirplot(sw.IRbar, opt, sw.IRinf, sw.IRsup)
    for (j, p) in enumerate(ps)
        out["SW2001_IR_$j"] = p
    end
    out["SW2001_VD"] = VARvdplot(sw.VDbar, opt)

    # ---- Blanchard & Quah (1989): VARirplot + hand-rolled Figure 1 ----
    _, _, bq_long, _ = load_replication_data(joinpath("BQ1989", "BQ1989_Data.xlsx"))
    bq = blanchard_quah_1989(ndraws = 1000,
                             band_src = DrawSource(uniforms = _pool(R["bq_irband_draws"])))
    opt = bq.VARopt
    opt = @set opt.vnames = bq_long
    opt = @set opt.FigSize = [26, 12]
    opt = @set opt.figname = joinpath(outdir, "BQ1989_")
    ps = VARirplot(bq.IRbar, opt, bq.IRinf, bq.IRsup)
    for (j, p) in enumerate(ps)
        out["BQ1989_IR_$j"] = p
    end

    # GO_BQ1989.m:65-88 — Figure 1 of the paper, from the raw (point) IR:
    # cumulated GDP + unemployment per shock, demand shock sign-flipped,
    # dashed zero on the left panel, solid on the right (as upstream).
    nsteps = size(bq.IR, 1)
    p = plot(layout = (1, 2), size = FigSize(26, 8))
    for sp in p.subplots
        E._matlab_axes!(sp)
    end
    plot!(p[1], 1:nsteps, cumsum(bq.IR[:, 1, 1]); linewidth = 2.5,
          linecolor = cmap(1), label = "GDP Level")
    plot!(p[1], 1:nsteps, bq.IR[:, 2, 1]; linewidth = 2.5,
          linecolor = cmap(2), label = "Unemployment")
    plot!(p[1], 1:nsteps, zeros(nsteps); linecolor = :black,
          linestyle = :dash, linewidth = 0.5, label = "")
    plot!(p[1]; title = "Supply shock", legend = :topright,
          legend_background_color = :white)
    plot!(p[2], 1:nsteps, cumsum(-bq.IR[:, 1, 2]); linewidth = 2.5,
          linecolor = cmap(1), label = "GDP Level")
    plot!(p[2], 1:nsteps, -bq.IR[:, 2, 2]; linewidth = 2.5,
          linecolor = cmap(2), label = "Unemployment")
    plot!(p[2], 1:nsteps, zeros(nsteps); linecolor = :black,
          linestyle = :solid, linewidth = 0.5, label = "")
    plot!(p[2]; title = "Demand shock", legend = :topright,
          legend_background_color = :white)
    SaveFigure(p, joinpath(outdir, "BQ_Replication"))
    out["BQ_Replication"] = p

    # ---- Uhlig (2005): hand-rolled PlotSwathe panels ------------------
    # GO_Uhlig2005.m:91-102 — 3×2 grid filled column-major (idx = [1 3 5 2
    # 4 6]), PlotSwathe with defaults on the median/68% bands of the single
    # restricted shock, titles from the SHORT variable names, axis tight.
    # The output filename's typo ("Uhligh") is upstream's, kept verbatim.
    _, uh_vnames, _, _ = load_replication_data(joinpath("Uhlig2005", "Uhlig2005_Data.xlsx"))
    uh = uhlig_2005(ndraws = Int(R["uh_ndraws"]),
                    src = DrawSource(normals = _pool(R["uh_sr_draws"])),
                    post_draws = _postpool(R["uh_post"]))
    SRout = uh.SRout
    nsteps = size(SRout.IRmed, 1)
    nvar = size(SRout.IRmed, 2)
    idx = [1, 3, 5, 2, 4, 6]
    p = plot(layout = (3, 2), size = FigSize(20, 24))
    for sp in p.subplots
        E._matlab_axes!(sp)
    end
    for ii in 1:nvar
        sp = p[idx[ii]]
        med = SRout.IRmed[:, ii, 1]
        inf = SRout.IRinf[:, ii, 1]
        sup = SRout.IRsup[:, ii, 1]
        PlotSwathe!(sp, med, hcat(inf, sup))
        plot!(sp, 1:nsteps, zeros(nsteps); linecolor = :black,
              linestyle = :dash, linewidth = 0.5, label = "")
        plot!(sp; title = uh_vnames[ii])
        _axis_tight!(sp, nsteps, med, inf, sup)
    end
    SaveFigure(p, joinpath(outdir, "Uhligh_Replication"))
    out["Uhligh_Replication"] = p

    # ---- Gertler & Karadi (2015): hand-rolled Cholesky-vs-IV grid -----
    # GO_GK2015.m:79-118 — 4×2: left column recursive (red, solid mean,
    # dashed bands), right column external-instrument (black), solid zero
    # line, axis tight, a "Cholesky"/"IV" legend on the bottom panel of
    # each column.
    gk = gertler_karadi_2015(ndraws = Int(R["gk_ndraws"]),
                             chol_src = DrawSource(uniforms = _pool(R["gk_chol_draws"])),
                             iv_src = DrawSource(uniforms = _pool(R["gk_iv_draws"])))
    gk_long = ["1yr rate", "CPI", "IP", "EBP"]
    nsteps = size(gk.cBAR, 1)
    p = plot(layout = (4, 2), size = FigSize(20, 20))
    for sp in p.subplots
        E._matlab_axes!(sp)
    end
    for (colidx, BAR, INF, SUP, colr, lab) in
        ((1, gk.cBAR, gk.cINF, gk.cSUP, :red, "Cholesky"),
         (2, gk.iBAR, gk.iINF, gk.iSUP, :black, "IV"))
        for ii in 1:4
            sp = p[colidx + 2 * (ii - 1)]
            islast = ii == 4
            plot!(sp, 1:nsteps, BAR[:, ii, 1]; linecolor = colr, linewidth = 2,
                  label = islast ? lab : "")
            plot!(sp, 1:nsteps, INF[:, ii, 1]; linecolor = colr, linewidth = 1,
                  linestyle = :dash, label = "")
            plot!(sp, 1:nsteps, SUP[:, ii, 1]; linecolor = colr, linewidth = 1,
                  linestyle = :dash, label = "")
            plot!(sp, 1:nsteps, zeros(nsteps); linecolor = :black,
                  linewidth = 0.5, label = "")
            plot!(sp; title = gk_long[ii])
            _axis_tight!(sp, nsteps, BAR[:, ii, 1], INF[:, ii, 1], SUP[:, ii, 1])
            islast && plot!(sp; legend = :topright, legend_background_color = :white)
        end
    end
    SaveFigure(p, joinpath(outdir, "GK_Replication"))
    out["GK_Replication"] = p

    # ---- Primer §9: SW historical decomposition (VARhdplot + dates) ---
    ps = primer_sw_section9(outdir)
    for (j, pp) in enumerate(ps)
        out["SW_HD_$j"] = pp
    end

    return out
end
