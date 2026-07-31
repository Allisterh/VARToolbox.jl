# =======================================================================
# Phase 6: end-to-end replication of the four published studies.
#
# The brief proposed comparing figures "by eye and by number". Comparison by
# eye is not evidence, so each replication is asserted elementwise against
# output dumped from the MATLAB original by gen_phase6.m -- including the
# bootstrap bands and the sign-restriction draws, which are exact rather
# than distributional thanks to the recorded draw sequences (M2).
# =======================================================================

include(joinpath(@__DIR__, "..", "examples", "replications.jl"))

R = fixture("phase6")
rupool(x) = [Matrix{Float64}(a) for a in vec(x)]

"MATLAB cell of {sigma_draw, Ft_draw} pairs -> Vector of tuples"
function postpool(x)
    [(Matrix{Float64}(p[1]), Matrix{Float64}(p[2])) for p in vec(x)]
end

@testset "Phase 6 -- data loading matches xlsread" begin
    # One data source for both languages: XLSX.jl must read the fixtures
    # identically to MATLAB's xlsread, sentinel handling included.
    X, _, _, _ = load_replication_data(joinpath("SW2001", "SW2001_Data.xlsx"))
    @golden X R["sw_X"]
    Xb, _, _, _ = load_replication_data(joinpath("BQ1989", "BQ1989_Data.xlsx"))
    @golden Xb R["bq_X"]
end

@testset "Phase 6 -- Stock & Watson (2001)" begin
    src = DrawSource(uniforms = rupool(R["sw_irband_draws"]))
    vsrc = DrawSource(uniforms = rupool(R["sw_vdband_draws"]))
    sw = stock_watson_2001(ndraws = Int(R["sw_ndraws"]),
                           band_src = src, vd_src = vsrc)

    @golden sw.IR R["sw_IR"]
    @golden sw.B  R["sw_B"]
    @golden sw.VD R["sw_VD"]

    @testset "bootstrap bands (exact, 1000 draws)" begin
        @golden sw.IRinf R["sw_IRinf"]
        @golden sw.IRsup R["sw_IRsup"]
        @golden sw.IRmed R["sw_IRmed"]
        @golden sw.IRbar R["sw_IRbar"]
        @golden sw.VDinf R["sw_VDinf"]
        @golden sw.VDsup R["sw_VDsup"]
        @golden sw.VDmed R["sw_VDmed"]
        @golden sw.VDbar R["sw_VDbar"]
    end

    @testset "Table 1.B" begin
        @golden sw.table1B R["sw_table1B"]
        # Each row is a variance decomposition and must sum to 100.
        for r in 1:size(sw.table1B, 1)
            @test sum(sw.table1B[r, :]) ≈ 100.0 rtol = 1e-9
        end
    end

    @testset "economic sanity" begin
        # Recursive ordering (unemp, infl, ff): a funds-rate shock cannot
        # move unemployment or inflation on impact.
        @test sw.IR[1, 1, 3] ≈ 0.0 atol = 1e-12
        @test sw.IR[1, 2, 3] ≈ 0.0 atol = 1e-12
        @test sw.B[1, 2] ≈ 0.0 atol = 1e-12
        @test sw.VAR.maxEig < 1.0
    end
end

@testset "Phase 6 -- Blanchard & Quah (1989)" begin
    src = DrawSource(uniforms = rupool(R["bq_irband_draws"]))
    bq = blanchard_quah_1989(ndraws = 1000, band_src = src)

    @golden bq.IR R["bq_IR"]
    @golden bq.B  R["bq_B"]
    @golden bq.IRinf R["bq_IRinf"]
    @golden bq.IRsup R["bq_IRsup"]
    @golden bq.IRmed R["bq_IRmed"]
    @golden bq.IRbar R["bq_IRbar"]

    @testset "Figure 2 series" begin
        @golden bq.fig2_supply_gdp R["bq_fig2_supply_gdp"]
        @golden bq.fig2_supply_u   R["bq_fig2_supply_u"]
        @golden bq.fig2_demand_gdp R["bq_fig2_demand_gdp"]
        @golden bq.fig2_demand_u   R["bq_fig2_demand_u"]
    end

    @testset "the long-run restriction actually binds" begin
        # This is the identifying assumption of the paper: the demand shock
        # has no permanent effect on the level of output. The cumulated GDP
        # response to shock 2 must die out.
        cum_demand = cumsum(bq.IR[:, 1, 2])
        @test abs(cum_demand[end]) < 0.02 * maximum(abs, cum_demand)
        # The supply shock, by contrast, has a permanent effect.
        cum_supply = cumsum(bq.IR[:, 1, 1])
        @test abs(cum_supply[end]) > 0.5 * maximum(abs, cum_supply)
    end
end

@testset "Phase 6 -- Gertler & Karadi (2015)" begin
    csrc = DrawSource(uniforms = rupool(R["gk_chol_draws"]))
    isrc = DrawSource(uniforms = rupool(R["gk_iv_draws"]))
    gk = gertler_karadi_2015(ndraws = Int(R["gk_ndraws"]),
                             chol_src = csrc, iv_src = isrc)

    # Conditioning tolerance, as in Phase 3: 49 coefficients on 384
    # persistent observations give cond(X'X) ~ 2.2e9, so the normal
    # equations cannot deliver 1e-10 in any implementation.
    IVTOL = 1e-6
    @test cond(gk.VAR.X' * gk.VAR.X) > 1e8

    @test approxeq(gk.IR_chol, R["gk_IR_chol"]; rtol = IVTOL)
    @test approxeq(gk.B_chol, R["gk_B_chol"]; rtol = IVTOL)
    @test approxeq(gk.cINF, R["gk_chol_INF"]; rtol = IVTOL)
    @test approxeq(gk.cSUP, R["gk_chol_SUP"]; rtol = IVTOL)
    @test approxeq(gk.cMED, R["gk_chol_MED"]; rtol = IVTOL)
    @test approxeq(gk.cBAR, R["gk_chol_BAR"]; rtol = IVTOL)

    @test approxeq(gk.IR_iv, R["gk_IR_iv"]; rtol = IVTOL)
    @test approxeq(gk.B_iv, R["gk_B_iv"]; rtol = IVTOL)
    @test approxeq(gk.Biv, R["gk_Biv"]; rtol = IVTOL)
    @test approxeq(gk.iINF, R["gk_iv_INF"]; rtol = IVTOL)
    @test approxeq(gk.iSUP, R["gk_iv_SUP"]; rtol = IVTOL)
    @test approxeq(gk.iMED, R["gk_iv_MED"]; rtol = IVTOL)
    @test approxeq(gk.iBAR, R["gk_iv_BAR"]; rtol = IVTOL)

    @testset "identification differs from the recursive benchmark" begin
        # If IV and Cholesky agreed, the instrument would be doing nothing.
        @test !approxeq(gk.IR_iv, gk.IR_chol; rtol = 1e-2)
        # iv fills only the first column of B.
        @test all(gk.B_iv[:, 2:end] .== 0)
    end
end

@testset "Phase 6 -- Uhlig (2005)" begin
    src = DrawSource(normals = rupool(R["uh_sr_draws"]))
    post = postpool(R["uh_post"])
    uh = uhlig_2005(ndraws = Int(R["uh_ndraws"]), src = src, post_draws = post)

    # sr_mod = 1, so each draw also redraws (sigma, F) from the posterior.
    # Those draws come from wishrnd/mvnrnd, whose streams cannot be
    # reproduced across languages -- so gen_phase6.m records them and they
    # are injected here. With both the rotation draws and the posterior
    # draws injected, this is an exact comparison.
    @test Int(R["uh_sr_mod"]) == 1

    @golden uh.SRout.Ball  R["uh_Ball"]
    @golden uh.SRout.B     R["uh_B"]
    @golden uh.SRout.Bmed  R["uh_Bmed"]
    @golden uh.SRout.IRmed R["uh_IRmed"]
    @golden uh.SRout.IRinf R["uh_IRinf"]
    @golden uh.SRout.IRsup R["uh_IRsup"]
    @golden uh.SRout.VDmed R["uh_VDmed"]

    @testset "every accepted draw satisfies the sign restrictions" begin
        SIGN = uh.SIGN
        for d in 1:size(uh.SRout.Ball, 3)
            B = uh.SRout.Ball[:, :, d]
            for i in 1:size(SIGN, 1)
                SIGN[i, 1] == 0 && continue
                @test sign(B[i, 1]) == sign(SIGN[i, 1])
            end
        end
    end

    @testset "the identified shock behaves like monetary policy" begin
        # Median impact responses must follow the imposed pattern: prices,
        # commodity prices and nonborrowed reserves down, funds rate up.
        m = uh.SRout.IRmed[1, :, 1]
        @test m[2] < 0    # deflator
        @test m[3] < 0    # commodity price
        @test m[5] < 0    # nonborrowed reserves
        @test m[6] > 0    # fed funds rate
    end

    @testset "bands bracket the median" begin
        @test all(uh.SRout.IRinf .<= uh.SRout.IRmed .+ 1e-12)
        @test all(uh.SRout.IRmed .<= uh.SRout.IRsup .+ 1e-12)
    end
end
