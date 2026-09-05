# =======================================================================
# Phase 1 golden-master tests: Utils, Auxiliary, and the OLS engine.
#
# Every comparison in this file is against output dumped from the MATLAB
# original by fixtures/gen/gen_phase1.m. No RNG is involved anywhere, so
# these are exact golden-master checks at rtol = 1e-10, not distributional
# comparisons.
# =======================================================================

F = fixture("phase1")
A = F["A"]
N = F["N"]

@testset "Phase 1 -- Utils" begin

    @testset "CommonSample" begin
        OUT, fo, lo = CommonSample(N)
        @golden OUT F["cs_out"]
        @test fo == Int(F["cs_fo"])

        # The fixture deliberately contains an INTERIOR NaN row. Under VAR
        # Toolbox 3.1 `lo` counted every NaN row outside the leading run
        # (the v3.1 fixture records 3); v4.0 fixed the loop so `lo` is the
        # trailing-run length. The port follows v4.0 -- see test_phase9.jl
        # for the v4.0 golden value. Pin the v3.1 value as NOT reproduced.
        @test Int(F["cs_lo"]) == 3
        @test lo == 2          # 2 trailing rows; the interior row no longer counts
        @test fo == 2

        OUT2, fo2, lo2 = CommonSample(permutedims(N), 2)
        @golden OUT2 F["cs2_out"]
        @test fo2 == Int(F["cs2_fo"])
        @test Int(F["cs2_lo"]) == 3 && lo2 == 2

        OUT3, fo3, lo3 = CommonSample(A)
        @golden OUT3 F["cs3_out"]
        @test (fo3, lo3) == (Int(F["cs3_fo"]), Int(F["cs3_lo"])) == (0, 0)
    end

    @testset "design matrices" begin
        @golden VARmakelags(A, 2) F["mkl_2"]
        @golden VARmakelags(A, 1) F["mkl_1"]
        for c in 0:3
            Y, X = VARmakexy(A, 2, c)
            @golden Y F["mkxy_Y_c$c"]
            @golden X F["mkxy_X_c$c"]
        end
    end

    @testset "lag / lead operators" begin
        @golden VARToolbox.lag(A, 1) F["lag_1"]
        @golden VARToolbox.lag(A, 2) F["lag_2"]
        @golden VARToolbox.lag(A, 2, -99) F["lag_2v"]
        @golden trimr(A, 2, 1) F["trimr_a"]
        @test_throws ArgumentError trimr(A, 5, 5)

        # L's positive argument is a LEAD -- upstream convention, preserved.
        Lp = L(A, 2)
        Ln = L(A, -2)
        @test isnan.(Lp[end-1:end, :]) == trues(2, 3)     # lead pads the tail
        @test isnan.(Ln[1:2, :]) == trues(2, 3)           # lag pads the head
        @test isequal(Lp, F["L_pos"])
        @test isequal(Ln, F["L_neg"])
        @golden L(A, 0) F["L_zero"]
    end

    @testset "misc" begin
        @golden seqa(2, 0.5, 7) F["seqa"]
        @golden ptrend(0, 8) F["ptrend_0"]
        @golden ptrend(2, 8) F["ptrend_2"]
        @golden getqr(A[1:3, 1:3]) F["getqr"]
        @golden getqr([2.0 -1 0; -1 2 -1; 0 -1 2]) F["getqr_b"]
        @test isequal(Num2NaN([123456789.0 2; 3 123456789; 5 6]), F["num2nan"])
        @golden NaN2Num([NaN 2.0; 3 NaN; 5 6]) F["nan2num"]

        # getqr's contract: orthonormal, with diag(R) >= 0 fixed by sign flips.
        Q = getqr(A[1:3, 1:3])
        @test Q * Q' ≈ I rtol = RTOL
    end
end

@testset "Phase 1 -- prctile (M1)" begin
    pv = flat(F["pv"])
    levels = flat(F["pct_levels"])

    @testset "matches MATLAB, not Julia's default quantile" begin
        @golden VARToolbox._prctile(pv, levels) F["pct_vec"]

        # The whole point of M1: Julia's default quantile does NOT match.
        # If this test ever starts passing, someone has changed the
        # definition and every error band in the package is now wrong.
        using Statistics: quantile
        naive = [quantile(pv, p / 100) for p in levels]
        @test !isapprox(naive, flat(F["pct_vec"]); rtol = 1e-6)
    end

    @testset "dimensions" begin
        @golden VARToolbox._prctile(A, [2.5, 50, 97.5], 1) F["pct_mat_d1"]
        @golden VARToolbox._prctile(A, [2.5, 50, 97.5], 2) F["pct_mat_d2"]
        # 4-D along dim 4 is what SR.m and VARirband.m actually do.
        R4 = F["R4"]
        got = VARToolbox._prctile(R4, [16, 50, 84], 4)
        @test size(got) == size(F["pct_4d"])
        @golden got F["pct_4d"]
    end

    @testset "NaN handling" begin
        @golden VARToolbox._prctile(flat(F["pn"]), levels) F["pct_nan"]
        @test isnan(VARToolbox._prctile([NaN, NaN], 50))
    end
end

@testset "Phase 1 -- distributions (M3)" begin
    tx = flat(F["tx"])
    tp = flat(F["tp"])

    @testset "tdis_prb" begin
        for n in (5, 30, 200)
            @golden tdis_prb(tx, n) F["tdis_prb_$n"]
        end

        # The upstream clamp (x2 -> 1-1e-12) perturbs the value at x=0 by
        # ~1e-6, because the Beta(n/2,1/2) density is singular at 1. The
        # closed form would give exactly 1. This asserts the clamp is
        # actually being replicated -- substituting 2*ccdf(TDist(n),|x|)
        # would break it.
        at_zero = tdis_prb([0.0], 5)[1]
        @test at_zero < 1.0
        @test 1.0 - at_zero > 1e-8
        @test 1.0 - at_zero < 1e-4
    end

    @testset "tdis_inv" begin
        for n in (5, 30, 200)
            @golden tdis_inv(tp, n) F["tdis_inv_$n"]
        end

        # Upstream loses ~7 digits at the centre of the distribution: the
        # final sqrt amplifies beta_inv's ~1e-14 Newton residual. The exact
        # answer here is 0; MATLAB returns 3.2148354647791181e-07. Pinned to
        # MATLAB's value, not the mathematically correct one -- if this ever
        # starts returning 0, the algorithm has been silently substituted.
        @test tdis_inv([0.5], 10)[1] ≈ 3.2148354647791181e-07 rtol = 1e-13
        @test tdis_inv([0.5], 10)[1] != 0.0
    end

    @testset "beta / F" begin
        @golden beta_pdf([0.01, 0.1, 0.5, 0.9, 0.99], 2, 3) F["beta_pdf_23"]
        @golden beta_inv([0.01, 0.1, 0.5, 0.9, 0.99], 2, 3) F["beta_inv_23"]
        @golden fdis_prb([0.1, 0.5, 1.0, 2.0, 5.0, 10.0], 3, 5) F["fdis_prb_35"]
        @test beta_pdf([-0.5, 1.5], 2, 3) == [0.0, 0.0]
    end
end

@testset "Phase 1 -- nan* family" begin
    M = F["M"]
    @golden nanmean(M, 1) F["nanmean_1"]
    @golden nanmean(M, 2) F["nanmean_2"]
    @golden nansum(M, 1) F["nansum_1"]
    @golden nansum(M, 2) F["nansum_2"]
    @golden nanstd(M, 1) F["nanstd_1"]
    @golden nanvar(M, 1) F["nanvar_1"]
    @golden nansem(M, 1) F["nansem_1"]

    # An all-NaN slice is NaN, not zero -- differs from sum(skipmissing(...)).
    @test isnan(nansum(M, 1)[3]) == isnan(F["nansum_1"][3])

    @testset "nanmedian zero-median defect" begin
        @golden nanmedian(M, 1) F["nanmedian_1"]
        # Upstream turns an exactly-zero median into NaN. Pinned so the
        # defect cannot be "fixed" without an intentional decision.
        Z = [-1.0 1.0; 1.0 3.0]
        @test isnan(nanmedian(Z, 1)[1])
        @test nanmedian(Z, 1)[2] == 2.0
    end
end

@testset "Phase 1 -- OLSmodel" begin
    y = A[:, 1]
    x = A[:, 2:3]

    for c in 0:2
        @testset "const_ = $c" begin
            O = OLSmodel(y, x, c)
            p = "ols_c$(c)_"
            @golden O.beta    F[p*"beta"]
            @golden O.yhat    F[p*"yhat"]
            @golden O.resid   F[p*"resid"]
            @golden O.bstd    F[p*"bstd"]
            @golden O.tstat   F[p*"tstat"]
            # tprob / bint: v3.1 used `nobs` degrees of freedom, v4.0 uses
            # `nobs - nvar`; the port follows v4.0 and these two fields are
            # checked against the v4.0 fixture in test_phase9.jl.
            @test !approxeq(O.tprob, F[p*"tprob"]; rtol = 1e-6)
            @golden O.sigbeta F[p*"sigbeta"]
            @golden O.bstd_HW F[p*"bstd_HW"]
            @golden O.bstd_NW F[p*"bstd_NW"]
            @test O.sige ≈ F[p*"sige"] rtol = RTOL
            @test O.rsqr ≈ F[p*"rsqr"] rtol = RTOL
            @test O.rbar ≈ F[p*"rbar"] rtol = RTOL
            @test O.dw   ≈ F[p*"dw"]   rtol = RTOL
            @test O.nvar == Int(F[p*"nvar"])
            @test O.nobs == Int(F[p*"nobs"])
            c > 0 && @test O.F ≈ F[p*"F"] rtol = RTOL
        end
    end

    @testset "const_ = 3 is rejected" begin
        # Upstream errors here (adds 2 columns, increments nvar by 3).
        @test F["ols_c3_errors"] == 1
        @test_throws ArgumentError OLSmodel(y, x, 3)
    end

    @testset "F is NaN when const_ == 0" begin
        @test isnan(OLSmodel(y, x, 0).F)
    end
end
