# =======================================================================
# Phase 5 golden-master tests: LPoption, LPmakex, LPmodel.
#
# Local projections involve no RNG, so every check here is an exact
# golden-master comparison at rtol = 1e-10.
# =======================================================================

Q = fixture("phase5")
SW5 = Q["SW"]

@testset "Phase 5 -- LPoption defaults" begin
    o = LPoption()
    @test o.nsteps == Int(Q["opt_nsteps"])
    @test o.pctg == Q["opt_pctg"]
    @test o.impact == Int(Q["opt_impact"])
    @test o.EXOG === nothing
    @test o.LAGS === nothing
    @test Q["opt_EXOG_empty"] == 1
    @test Q["opt_LAGS_empty"] == 1

    @testset "LPoption.m omits nlag and const" begin
        # Pinned as a contract: the MATLAB options struct genuinely lacks
        # these fields, so LPmodel.m dies on a non-existent-field access.
        @test Q["opt_has_nlag"] == 0
        @test Q["opt_has_const"] == 0
        @test Q["lp_no_nlag_errors"] == 1
        # This port declares them, with nlag unset by default and an
        # informative error rather than a field-access failure.
        @test o.nlag === nothing
        @test o.const_ == 1
        opt = @set o.LAGS = SW5
        @test_throws ArgumentError LPmodel(SW5[:, 1:1], SW5[:, 3:3], opt)
    end
end

@testset "Phase 5 -- LPmakex" begin
    A = Q["A"]
    @golden LPmakex(A, 1) Q["lpmakex_1"]
    @golden LPmakex(A, 3) Q["lpmakex_3"]
    @golden LPmakex(A[:, 1], 2) Q["lpmakex_1col"]

    @testset "lags, not leads" begin
        # LPmakex builds on L with a NEGATIVE argument, because in this
        # toolbox a positive argument to L is a lead.
        X = LPmakex(A, 2)
        @test all(isnan, X[1, 1:size(A, 2)])           # first lag: 1 NaN row
        @test all(isnan, X[1:2, (size(A, 2)+1):end])   # second lag: 2 NaN rows
        @test X[2, 1:size(A, 2)] ≈ A[1, :] rtol = RTOL
    end
end

ENDO5 = SW5[:, 1:1]
S5 = SW5[:, 3:3]

@testset "Phase 5 -- LPmodel with LAGS" begin
    opt = LPoption()
    opt = @set opt.nsteps = Int(Q["lp1_nsteps"])
    opt = @set opt.nlag = Int(Q["lp1_nlag"])
    opt = @set opt.const_ = 1
    opt = @set opt.LAGS = SW5
    opt = @set opt.impact = 0

    IR, INF, SUP, LPout = LPmodel(ENDO5, S5, opt)
    @golden IR  Q["lp1_IR"]
    @golden INF Q["lp1_INF"]
    @golden SUP Q["lp1_SUP"]
    @test LPout.nvar == Int(Q["lp1_nvar"])
    @test LPout.ntotcoeff == Int(Q["lp1_ntotcoeff"])
    @test length(LPout.h) == opt.nsteps

    for hh in 1:opt.nsteps
        q = "lp1_h$(hh)_"
        a = LPout.h[hh]
        @test a.nobs == Int(Q[q*"nobs"])
        @golden a.beta    Q[q*"beta"]
        @golden a.tstat   Q[q*"tstat"]
        @golden a.bstd    Q[q*"bstd"]
        @golden a.bstd_HW Q[q*"bstd_HW"]
        @golden a.tprob   Q[q*"tprob"]
        @test a.rsqr  ≈ Q[q*"rsqr"]  rtol = RTOL
        @test a.dw    ≈ Q[q*"dw"]    rtol = RTOL
        @test a.sigma ≈ Q[q*"sigma"] rtol = RTOL
    end

    @testset "IR[h] is the shock coefficient" begin
        for hh in 1:opt.nsteps
            @test IR[hh] == LPout.h[hh].beta[1]
        end
    end

    @testset "bstd_HW actually holds Newey-West values" begin
        # LPmodel.m assigns bstd_HW twice in a row, the second from
        # OLSout.bstd_NW, so the stored field is mislabelled. Both are
        # exposed here; this pins that the MATLAB field matches NW and not
        # White, and that the two genuinely differ.
        a = LPout.h[1]
        X = a.x
        O = OLSmodel(a.y, X, 0)
        @test approxeq(a.bstd_HW, O.bstd_NW)
        @test approxeq(a.bstd_White, O.bstd_HW)
        @test !approxeq(O.bstd_NW, O.bstd_HW; rtol = 1e-6)
    end

    @testset "bands use the WHITE errors, not the stored field" begin
        # INF/SUP come from OLSout.bstd_HW (White), while LPout.h[].bstd_HW
        # stores Newey-West. So the bands cannot be reconstructed from the
        # stored field -- a genuine inconsistency, pinned.
        conf = quantile(Normal(), 1 - (1 - opt.pctg / 100) / 2)
        from_white = IR .- conf .* [LPout.h[hh].bstd_White[1] for hh in 1:opt.nsteps]
        @test approxeq(from_white, INF)
        from_nw = IR .- conf .* [LPout.h[hh].bstd_HW[1] for hh in 1:opt.nsteps]
        @test !approxeq(from_nw, INF; rtol = 1e-6)
    end

    @testset "bands are symmetric about IR" begin
        @test approxeq(SUP .- IR, IR .- INF)
    end
end

@testset "Phase 5 -- LPmodel unit shock" begin
    opt = LPoption()
    opt = @set opt.nsteps = Int(Q["lp1_nsteps"])
    opt = @set opt.nlag = Int(Q["lp1_nlag"])
    opt = @set opt.LAGS = SW5
    opt = @set opt.impact = 1
    IR, INF, SUP, _ = LPmodel(ENDO5, S5, opt)
    @golden IR  Q["lp2_IR"]
    @golden INF Q["lp2_INF"]
    @golden SUP Q["lp2_SUP"]

    @testset "impact=0 rescales impact=1 by the shock's std" begin
        o0 = @set opt.impact = 0
        IR0, _, _, _ = LPmodel(ENDO5, S5, o0)
        sd = std(S5[(1+opt.nlag):end, 1])
        @test approxeq(IR0, IR .* sd; rtol = 1e-9)
    end
end

@testset "Phase 5 -- LPmodel without LAGS" begin
    opt = LPoption()
    opt = @set opt.nsteps = Int(Q["lp3_nsteps"])
    opt = @set opt.nlag = 4
    opt = @set opt.const_ = 1
    opt = @set opt.LAGS = nothing
    IR, INF, SUP, LPout = LPmodel(ENDO5, S5, opt)
    @golden IR  Q["lp3_IR"]
    @golden INF Q["lp3_INF"]
    @golden SUP Q["lp3_SUP"]
    @test LPout.ntotcoeff == Int(Q["lp3_ntotcoeff"])
    @golden LPout.h[1].beta Q["lp3_h1_beta"]
    @test LPout.h[1].nobs == Int(Q["lp3_h1_nobs"])

    @testset "no intercept when LAGS is absent" begin
        # The constant enters only through VARmakexy on LAGS. With no LAGS
        # the projection has no intercept at all, whatever const_ says.
        @test size(LPout.h[1].x, 2) == 1        # the shock, and nothing else
        # ...and const_ genuinely does not change the fit.
        o2 = @set opt.const_ = 0
        IR2, _, _, _ = LPmodel(ENDO5, S5, o2)
        @test approxeq(IR2, IR)
    end
end

@testset "Phase 5 -- LPmodel with EXOG" begin
    EX = reshape(flat(Q["EX"]), :, 1)
    opt = LPoption()
    opt = @set opt.nsteps = Int(Q["lp4_nsteps"])
    opt = @set opt.nlag = Int(Q["lp1_nlag"])
    opt = @set opt.LAGS = SW5
    opt = @set opt.EXOG = EX
    IR, INF, SUP, LPout = LPmodel(ENDO5, S5, opt)
    @golden IR  Q["lp4_IR"]
    @golden INF Q["lp4_INF"]
    @golden SUP Q["lp4_SUP"]
    @test LPout.ntotcoeff == Int(Q["lp4_ntotcoeff"])
    @golden LPout.h[1].beta Q["lp4_h1_beta"]
end

@testset "Phase 5 -- LPmodel coverage level" begin
    opt = LPoption()
    opt = @set opt.nsteps = Int(Q["lp5_nsteps"])
    opt = @set opt.nlag = Int(Q["lp1_nlag"])
    opt = @set opt.LAGS = SW5
    opt = @set opt.pctg = 68.0
    IR, INF, SUP, _ = LPmodel(ENDO5, S5, opt)
    @golden IR  Q["lp5_IR"]
    @golden INF Q["lp5_INF"]
    @golden SUP Q["lp5_SUP"]

    # A 68% band sits strictly inside the 95% band.
    o95 = @set opt.pctg = 95.0
    _, I95, S95, _ = LPmodel(ENDO5, S5, o95)
    @test all(INF .> I95)
    @test all(SUP .< S95)
end

@testset "Phase 5 -- error contracts" begin
    opt = LPoption()
    opt = @set opt.nsteps = 4
    opt = @set opt.nlag = 4
    opt = @set opt.LAGS = SW5

    @test Q["lp_bad_S_errors"] == 1
    @test_throws DimensionMismatch LPmodel(ENDO5, S5[1:(end-5), :], opt)

    @test Q["lp_bad_impact_errors"] == 1
    @test_throws ArgumentError LPmodel(ENDO5, S5, @set opt.impact = 7)

    @test_throws ArgumentError LPmodel(ENDO5, S5, @set opt.nlag = 0)
end
