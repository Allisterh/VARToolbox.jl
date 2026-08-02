# =======================================================================
# Phase 8: the plotting extension (ext/VARToolboxPlotsExt).
#
# Rendering cannot be pixel-asserted across GR versions and platforms, so
# these tests pin what IS assertable headlessly:
#   - the extension loads and every stub has a method;
#   - the primitives reproduce the Figure/ layer's constants and layout
#     logic (cmap colours, DatesPlot tick positions vs MATLAB, stacking);
#   - each plot function writes the upstream-named PDF files and its Plot
#     object carries exactly the input data (median = IR, band = INF/SUP);
#   - the three inconsistent `pick` meanings behave as upstream;
#   - examples/figures.jl regenerates all 12 replication figures.
#
# The visual comparison against MATLAB-rendered references is a separate,
# manual gate recorded in PORTING_LOG.md ("Phase 8").
# =======================================================================

ENV["GKSwstype"] = "100"    # headless GR
using Plots

@testset "Phase 8 -- plotting extension" begin

EXT = Base.get_extension(VARToolbox, :VARToolboxPlotsExt)

@testset "extension loads under the test environment" begin
    @test EXT !== nothing
    for f in (VARirplot, VARvdplot, VARhdplot, SRirplot, SRvdplot, SRhdplot,
              PlotSwathe, PlotSwathe!, AreaPlot!, cmap, DatesPlot!,
              FigSize, SaveFigure)
        @test length(methods(f)) >= 1
    end
end

@testset "cmap matches Figure/cmap.m exactly" begin
    C = [0.000 0.447 0.741
         0.929 0.694 0.125
         0.494 0.184 0.556
         0.301 0.745 0.933
         0.850 0.325 0.098
         0.466 0.674 0.188
         0.635 0.078 0.184
         0/255 48/255 143/255
         209/255 96/255 1/255]
    for i in 1:9
        c = cmap(i)
        @test (c.r, c.g, c.b) == (C[i, 1], C[i, 2], C[i, 3])
    end
    # The parula tail is not ported (PORTING_LOG Phase 8).
    @test_throws ErrorException cmap(10)
end

@testset "FigSize: cm -> px == PDF points" begin
    # GR maps `size` pixels 1:1 to PDF points (probed; MediaBox == size),
    # so px = cm * 72 / 2.54 gives a page physically matching MATLAB's
    # cm-sized figure window.
    @test FigSize(26, 12) == (737, 340)
    @test FigSize() == (737, 680)          # MATLAB FigSize dflt [26, 24]
    @test FigSize(20, 24) == (567, 680)
end

@testset "_date_ticks against MATLAB DatesPlot" begin
    # Expected values dumped from MATLAB R2026a (see PORTING_LOG Phase 8):
    #   figure; plot(1:nobs); DatesPlot(fo, nobs, nticks, freq);
    #   get(gca,'XTick'), get(gca,'XTickLabel'), get(gca,'XLim')
    # Case 1 -- the Primer §9 HD axis: quarterly from 1960Q1, 166 obs.
    pos, labs, lims = EXT._date_ticks(1960.0, 166, 8, "q")
    @test pos == [4.0, 24.0, 44.0, 64.0, 84.0, 104.0, 124.0, 144.0]
    @test labs == ["1960", "1965", "1970", "1975", "1980", "1985",
                   "1990", "1995"]
    @test lims == (0.0, 167.0)
    # Case 2 -- quarterly from a Q4 start (fo_diff = 1); upstream places
    # ticks past the data (25, 29) and lets XLim clip them.
    pos2, labs2, lims2 = EXT._date_ticks(1959.75, 24, 8, "q")
    @test pos2 == [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0]
    @test labs2 == ["1959", "1960", "1961", "1962", "1963", "1964",
                    "1965", "1966"]
    @test lims2 == (0.0, 25.0)
    # Case 3 -- annual.
    pos3, labs3, lims3 = EXT._date_ticks(2000, 20, 5, "y")
    @test pos3 == [4.0, 8.0, 12.0, 16.0, 20.0]
    @test labs3 == ["2003", "2007", "2011", "2015", "2019"]
    @test lims3 == (0.0, 21.0)
    # Case 4 -- monthly. Upstream's literal 0.0833 step (not 1/12) makes
    # fo_diff -- and so every tick position -- fractional: MATLAB places
    # ticks at 12.0048, 24.0048, ... and floors them when indexing the
    # label vector (with a console warning). Both are replicated.
    pos4, labs4, _ = EXT._date_ticks(1965.0, 100, 8, "m")
    @test length(pos4) == 7
    @test pos4 ≈ [12.0048, 24.0048, 36.0048, 48.0048, 60.0048, 72.0048,
                  84.0048] atol = 1e-3
    @test labs4 == ["1965", "1966", "1967", "1968", "1969", "1970", "1971"]
    @test_throws ErrorException EXT._date_ticks(2000.0, 20, 5, "w")
end

@testset "PlotSwathe series contract" begin
    bar = [1.0, 2.0, 3.0, 2.0]
    inf_ = bar .- 1.0
    sup_ = bar .+ 1.0
    p = PlotSwathe(bar, hcat(inf_, sup_))
    s = p[1].series_list
    # MATLAB net z-order: patch, low edge, upp edge, median on top. NB
    # PlotSwathe.m:85-89 names the INF column "upp" and the SUP column
    # "low" -- an upstream swap that is harmless for a filled band and is
    # replicated: the fill series' y is SUP with fillrange INF.
    @test length(s) == 4
    @test vec(s[1][:y]) == sup_
    @test vec(s[1][:fillrange]) == inf_
    @test s[1][:fillalpha] == 1.0            # opaque, as upstream
    @test s[1][:fillcolor] == Plots.RGBA(138 / 255, 178 / 255, 212 / 255, 1.0)
    @test vec(s[2][:y]) == sup_
    @test vec(s[3][:y]) == inf_
    @test vec(s[4][:y]) == bar
    @test s[4][:linecolor] == Plots.RGBA(13 / 255, 54 / 255, 84 / 255, 1.0)

    # Vector swathe means a symmetric band bar ± swathe.
    p2 = PlotSwathe(bar, fill(0.5, 4))
    s2 = p2[1].series_list
    @test vec(s2[1][:y]) == bar .- 0.5
    @test vec(s2[1][:fillrange]) == bar .+ 0.5

    # NaN band points are stripped from the patch, as MATLAB's does.
    inf_nan = copy(inf_); inf_nan[2] = NaN
    p3 = PlotSwathe(bar, hcat(inf_nan, sup_))
    @test length(vec(p3[1].series_list[1][:y])) == 3
end

@testset "AreaPlot! signed stacking" begin
    X = [1.0 2.0; -1.0 1.0; 2.0 -3.0]
    p = plot()
    AreaPlot!(p, X; labels = ["a", "b"])
    s = p[1].series_list
    @test length(s) == 4                     # pos stack then neg stack
    @test vec(s[1][:y]) == [1.0, 0.0, 2.0]   # cumsum of positive part, k=1
    @test vec(s[2][:y]) == [3.0, 1.0, 2.0]   # k=2 stacked on k=1
    @test vec(s[3][:y]) == [0.0, -1.0, 0.0]  # negative part, k=1
    @test vec(s[4][:y]) == [0.0, -1.0, -3.0]
    @test s[1][:fillcolor] == Plots.RGBA(cmap(1))
    @test s[2][:fillcolor] == Plots.RGBA(cmap(2))
    @test s[1][:label] == "a"                # positive stack carries labels
    @test s[3][:label] == ""
end

# ---- the six plot functions on a small, fixed VAR ---------------------

_pdf_ok(path) = isfile(path) && filesize(path) > 1000 &&
                read(path, 4) == UInt8['%', 'P', 'D', 'F']

rng8 = MersenneTwister(2026)
Y8 = zeros(140, 3)
for t in 2:140
    Y8[t, :] = 0.6 .* Y8[t - 1, :] .+ randn(rng8, 3)
end
VAR8, _ = VARmodel(Y8, 2)
opt8 = VARoption()
opt8 = @set opt8.nsteps = 12
opt8 = @set opt8.vnames = ["A", "B", "C"]
IR8, VAR8 = VARir(VAR8, opt8)
VD8, VAR8 = VARvd(VAR8, opt8)
HD8, VAR8 = VARhd(VAR8, opt8)
INF8 = IR8 .- 0.1
SUP8 = IR8 .+ 0.1

@testset "VARirplot: files, layout, data" begin
    dir = mktempdir()
    opt = @set opt8.figname = joinpath(dir, "t_")
    ps = VARirplot(IR8, opt, INF8, SUP8)
    @test length(ps) == 3                    # one figure per shock
    for j in 1:3
        @test _pdf_ok(joinpath(dir, "t_IR_$j.pdf"))
    end
    p = ps[1]
    @test length(p.subplots) == 4            # 2x2 grid, 4th blank
    @test isempty(p[4].series_list)
    @test p[1][:title] == "A to A"
    @test p[2][:title] == "B to A"
    for ii in 1:3
        s = p[ii].series_list
        # black under-line, patch, 2 edges, navy median, zero line.
        # The patch carries SUP as y and INF as fillrange -- the replicated
        # upstream naming swap, see the PlotSwathe testset.
        @test length(s) == 6
        @test vec(s[1][:y]) == IR8[:, ii, 1]
        @test vec(s[2][:y]) == SUP8[:, ii, 1]
        @test vec(s[2][:fillrange]) == INF8[:, ii, 1]
        @test vec(s[5][:y]) == IR8[:, ii, 1]
        @test all(iszero, s[6][:y])
    end
    # Without bands: no swathe series.
    dir2 = mktempdir()
    opt2 = @set opt8.figname = joinpath(dir2, "n_")
    ps2 = VARirplot(IR8, opt2)
    @test length(ps2[1][1].series_list) == 2
    # vnames are required, as upstream.
    @test_throws ErrorException VARirplot(IR8, VARoption())
end

@testset "VARvdplot: single file, area panels, legend row" begin
    dir = mktempdir()
    opt = @set opt8.figname = joinpath(dir, "t_")
    p = VARvdplot(VD8, opt)
    @test _pdf_ok(joinpath(dir, "t_VD.pdf"))
    # 2x2 grid + the bottom legend row (quality dflt 2).
    @test length(p.subplots) == 5
    s = p[1].series_list
    @test length(s) == 6                     # 3 pos + 3 neg stack series
    @test vec(s[1][:y]) ≈ VD8[:, 1, 1]       # k=1 stack == first shock share
    @test vec(s[3][:y]) ≈ vec(sum(VD8[:, :, 1], dims = 2))  # top of stack = 100
end

@testset "VARhdplot: per-variable files, data line, dates" begin
    dir = mktempdir()
    opt = @set opt8.figname = joinpath(dir, "t_")
    opt = @set opt.firstdate = 1980.0
    ps = VARhdplot(HD8, opt)
    @test length(ps) == 3
    for j in 1:3
        @test _pdf_ok(joinpath(dir, "t_HD_$j.pdf"))
    end
    p = ps[1]
    # area stacks + the black "Data" line == row sums of HD.shock
    s = p[1].series_list
    dataline = vec(s[7][:y])
    want = vec(sum(HD8.shock[:, :, 1], dims = 2))
    @test isequal(dataline, want)            # NaN-padded rows compare equal
    # DatesPlot applied: ticks are set, not :auto
    @test p[1][:xaxis][:ticks] != :auto
end

@testset "pick: three inconsistent upstream meanings" begin
    # VARirplot: pick = k -> ONLY shock k.
    dir = mktempdir()
    opt = @set opt8.figname = joinpath(dir, "t_")
    opt = @set opt.pick = 2
    ps = VARirplot(IR8, opt, INF8, SUP8)
    @test length(ps) == 1
    @test _pdf_ok(joinpath(dir, "t_IR_2.pdf"))
    @test !isfile(joinpath(dir, "t_IR_1.pdf"))
    # VARvdplot: pick = k -> variables 1..k.
    p = VARvdplot(VD8, opt)
    @test !isempty(p[1].series_list)
    @test !isempty(p[2].series_list)
    @test isempty(p[3].series_list)
    # VARhdplot: pick = k -> variables k..nvar.
    ps = VARhdplot(HD8, opt)
    @test length(ps) == 2
    @test _pdf_ok(joinpath(dir, "t_HD_2.pdf"))
    @test _pdf_ok(joinpath(dir, "t_HD_3.pdf"))
    @test !isfile(joinpath(dir, "t_HD_1.pdf"))
    # Out-of-range pick errors, as upstream.
    optbad = @set opt8.pick = 5
    @test_throws ErrorException VARirplot(IR8, optbad, INF8, SUP8)
end

@testset "SR variants require snames" begin
    dir = mktempdir()
    opt = @set opt8.figname = joinpath(dir, "t_")
    @test_throws ErrorException SRirplot(IR8, opt, INF8, SUP8)
    # nshocks comes from snames, not size(IR, 3): one name -> one figure.
    opt = @set opt.snames = ["Shock 1"]
    ps = SRirplot(IR8, opt, INF8, SUP8)
    @test length(ps) == 1
    @test _pdf_ok(joinpath(dir, "t_IR_1.pdf"))
    p = SRvdplot(VD8, opt)
    @test _pdf_ok(joinpath(dir, "t_VD.pdf"))
    ps = SRhdplot(HD8, opt)
    @test length(ps) == 3
end

@testset "examples/figures.jl regenerates all 12 figures" begin
    include(joinpath(@__DIR__, "..", "examples", "figures.jl"))
    dir = mktempdir()
    out = generate_all_figures(dir)
    stems = ["SW2001_IR_1", "SW2001_IR_2", "SW2001_IR_3", "SW2001_VD",
             "BQ1989_IR_1", "BQ1989_IR_2", "BQ_Replication",
             "Uhligh_Replication", "GK_Replication",
             "SW_HD_1", "SW_HD_2", "SW_HD_3"]
    @test sort(collect(keys(out))) == sort(stems)
    for s in stems
        @test _pdf_ok(joinpath(dir, "$s.pdf"))
    end
end

end # Phase 8
