# =======================================================================
# Phase 4 golden-master tests: VARirband, VARvdband.
#
# EXACT checks, not distributional. gen_phase4.m records the rand()
# sequence MATLAB consumed and this file replays it through a DrawSource
# (M2), so the bootstrap bands are compared elementwise at rtol = 1e-10.
#
# The brief originally proposed comparing summary statistics across 10,000
# draws at 1-2% relative tolerance. That would have been looser than the
# prctile defect M1 identifies, and would not have caught it.
# =======================================================================

P = fixture("phase4")
upool(x) = [Matrix{Float64}(a) for a in vec(x)]

SW4 = P["SW"]
V4, o4 = VARmodel(SW4, 4, 1)
o4 = @set o4.nsteps = Int(P["nsteps"])
o4 = @set o4.ndraws = Int(P["ndraws"])
o4 = @set o4.pctg = P["pctg"]

@testset "Phase 4 -- VARirband" begin
    @testset "standard bootstrap" begin
        src = DrawSource(uniforms = upool(P["irb_bs_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "bs"
        INF, SUP, MED, BAR = VARirband(V4, o; src = src)
        @golden INF P["irb_bs_INF"]
        @golden SUP P["irb_bs_SUP"]
        @golden MED P["irb_bs_MED"]
        @golden BAR P["irb_bs_BAR"]
        @test size(INF) == (Int(P["nsteps"]), V4.nvar, V4.nvar)
    end

    @testset "wild bootstrap" begin
        src = DrawSource(uniforms = upool(P["irb_wild_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "wild"
        INF, SUP, MED, BAR = VARirband(V4, o; src = src)
        @golden INF P["irb_wild_INF"]
        @golden SUP P["irb_wild_SUP"]
        @golden MED P["irb_wild_MED"]
        @golden BAR P["irb_wild_BAR"]
    end

    @testset "long-run identification" begin
        src = DrawSource(uniforms = upool(P["irb_long_draws"]))
        o = @set o4.ident = "long"
        o = @set o.method = "bs"
        INF, SUP, MED, BAR = VARirband(V4, o; src = src)
        @golden INF P["irb_long_INF"]
        @golden SUP P["irb_long_SUP"]
        @golden MED P["irb_long_MED"]
        @golden BAR P["irb_long_BAR"]
    end

    @testset "pctg = 68" begin
        src = DrawSource(uniforms = upool(P["irb_p68_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "bs"
        o = @set o.pctg = 68.0
        INF, SUP, _, _ = VARirband(V4, o; src = src)
        @golden INF P["irb_p68_INF"]
        @golden SUP P["irb_p68_SUP"]
        # A 68% band must sit inside the 95% band from the same draws.
        src95 = DrawSource(uniforms = upool(P["irb_p68_draws"]))
        o95 = @set o.pctg = 95.0
        I95, S95, _, _ = VARirband(V4, o95; src = src95)
        @test all(INF .>= I95 .- 1e-12)
        @test all(SUP .<= S95 .+ 1e-12)
    end

    @testset "with exogenous variable" begin
        EX = reshape(flat(P["EX"]), :, 1)
        Vx, ox = VARmodel(SW4, 4, 1, EX, 0)
        ox = @set ox.nsteps = Int(P["exog_nsteps"])
        ox = @set ox.ndraws = Int(P["exog_ndraws"])
        ox = @set ox.ident = "short"
        ox = @set ox.method = "bs"
        src = DrawSource(uniforms = upool(P["irb_exog_draws"]))
        INF, SUP, MED, BAR = VARirband(Vx, ox; src = src)
        @golden INF P["irb_exog_INF"]
        @golden SUP P["irb_exog_SUP"]
        @golden MED P["irb_exog_MED"]
        @golden BAR P["irb_exog_BAR"]
    end

    @testset "bands bracket the median" begin
        src = DrawSource(uniforms = upool(P["irb_bs_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "bs"
        INF, SUP, MED, _ = VARirband(V4, o; src = src)
        @test all(INF .<= MED .+ 1e-12)
        @test all(MED .<= SUP .+ 1e-12)
    end

    @testset "unknown method raises" begin
        @test P["band_bad_method_errors"] == 1
        o = @set o4.method = "nonsense"
        @test_throws ArgumentError VARirband(V4, o)
    end
end

@testset "Phase 4 -- VARirband with iv + wild" begin
    # This is the branch VARvdband does not have: the wild bootstrap
    # resamples the instrument alongside the residuals.
    ENDO = P["gk_ENDO"]
    IVv = P["gk_IV"]
    Vg, og = VARmodel(ENDO, 12, 1)
    Vg = @set Vg.IV = IVv
    og = @set og.nsteps = Int(P["iv_nsteps"])
    og = @set og.ndraws = Int(P["iv_ndraws"])
    og = @set og.ident = "iv"
    og = @set og.method = "wild"

    src = DrawSource(uniforms = upool(P["irb_iv_draws"]))
    INF, SUP, MED, BAR = VARirband(Vg, og; src = src)

    # Same conditioning justification as the Phase 3 IV block: cond(X'X) is
    # about 2.2e9 for this specification, so 1e-10 is unreachable through
    # the normal equations. See test_phase3.jl for the full argument.
    IVTOL = 1e-6
    @test approxeq(INF, P["irb_iv_INF"]; rtol = IVTOL)
    @test approxeq(SUP, P["irb_iv_SUP"]; rtol = IVTOL)
    @test approxeq(MED, P["irb_iv_MED"]; rtol = IVTOL)
    @test approxeq(BAR, P["irb_iv_BAR"]; rtol = IVTOL)
end

@testset "Phase 4 -- VARvdband" begin
    @testset "standard bootstrap" begin
        src = DrawSource(uniforms = upool(P["vdb_bs_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "bs"
        INF, SUP, MED, BAR = VARvdband(V4, o; src = src)
        @golden INF P["vdb_bs_INF"]
        @golden SUP P["vdb_bs_SUP"]
        @golden MED P["vdb_bs_MED"]
        @golden BAR P["vdb_bs_BAR"]
    end

    @testset "wild bootstrap" begin
        src = DrawSource(uniforms = upool(P["vdb_wild_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "wild"
        INF, SUP, MED, BAR = VARvdband(V4, o; src = src)
        @golden INF P["vdb_wild_INF"]
        @golden SUP P["vdb_wild_SUP"]
        @golden MED P["vdb_wild_MED"]
        @golden BAR P["vdb_wild_BAR"]
    end

    @testset "variance shares stay in [0, 100]" begin
        src = DrawSource(uniforms = upool(P["vdb_bs_draws"]))
        o = @set o4.ident = "short"
        o = @set o.method = "bs"
        INF, SUP, MED, BAR = VARvdband(V4, o; src = src)
        for A in (INF, SUP, MED, BAR)
            @test all(A .>= -1e-9)
            @test all(A .<= 100 + 1e-9)
        end
        # BAR is a mean across draws, each of which sums to 100 by shock.
        for h in 1:Int(P["nsteps"]), k in 1:V4.nvar
            @test sum(BAR[h, :, k]) ≈ 100.0 rtol = 1e-9
        end
    end
end

@testset "Phase 4 -- DrawSource contract" begin
    # An exhausted injected pool must raise, not silently fall back to
    # sampling -- a golden-master test that quietly stopped being one would
    # be worse than no test at all.
    src = DrawSource(uniforms = [rand(4, 1)])
    o = @set o4.ident = "short"
    o = @set o.method = "bs"
    o = @set o.ndraws = 5
    @test_throws ErrorException VARirband(V4, o; src = src)

    # Normal and uniform pools are tracked independently, so a pool
    # recorded from randn is never consumed by a rand call site: with only
    # `normals` injected, uniform draws still come from the rng and the
    # normal pool is left untouched.
    s2 = DrawSource(rng = MersenneTwister(9), normals = [randn(3, 3)])
    u = VARToolbox.next_uniform!(s2, 3, 1)
    @test size(u) == (3, 1)
    @test all(0 .<= u .< 1)
    @test s2.i == 0          # normal pool untouched
    @test s2.j == 0          # and no uniform pool to advance
    @test size(VARToolbox.next_normal!(s2, 3, 3)) == (3, 3)
    @test s2.i == 1          # normal pool advanced only by a normal draw
end
