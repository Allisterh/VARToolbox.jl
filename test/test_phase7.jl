# =======================================================================
# Stats/ golden-master tests. Deterministic throughout; rtol = 1e-10.
# =======================================================================

T7 = fixture("phase7")
X7 = T7["X"]
XN7 = T7["XN"]

@testset "Stats -- moving-window transforms" begin
    @golden MovAvg(X7, 5)  T7["movavg_5"]
    @golden MovAvg(X7, 12) T7["movavg_12"]
    @golden MovSum(X7, 5)  T7["movsum_5"]
    @golden MovStd(X7, 5)  T7["movstd_5"]
    @golden MovAvgCent(X7, 4) T7["movavgcent_4"]

    @testset "NaN handling differs between trailing and centred" begin
        # MovAvg uses a plain mean (NaN propagates); MovAvgCent uses
        # nanmean (NaN skipped). Two otherwise-parallel functions, two
        # conventions. Pinned so the inconsistency is visible.
        @golden MovAvg(XN7, 5)      T7["movavg_5_nan"]
        @golden MovAvgCent(XN7, 4)  T7["movavgcent_4_nan"]

        a = MovAvg(XN7, 5)
        c = MovAvgCent(XN7, 4)
        @test any(isnan, a[6:end, 2])          # propagated
        @test !any(isnan, c[10:60, 2])         # skipped
    end

    @testset "leading NaN padding" begin
        w = 5
        m = MovAvg(X7, w)
        @test all(isnan, m[1:(w-1), :])
        @test !any(isnan, m[w:end, :])
        @test size(m) == size(X7)
        # Row w is the mean of rows 1:w.
        @test m[w, 1] ≈ sum(X7[1:w, 1]) / w rtol = RTOL
    end

    @testset "window is a half-width in the centred versions" begin
        # MovAvgCent(X, w) averages 2w+1 rows, while MovAvg(X, w) averages
        # w. Upstream inconsistency, preserved.
        w = 4
        c = MovAvgCent(X7, w)
        @test all(isnan, c[1:w, :])
        @test all(isnan, c[(end-w+1):end, :])
        @test c[w+1, 1] ≈ mean(X7[1:(2w+1), 1]) rtol = RTOL
    end

    @testset "moving correlations" begin
        @golden MovCorr(X7, 10, 1)     T7["movcorr_10_1"]
        @golden MovCorr(X7, 10, 3)     T7["movcorr_10_3"]
        @golden MovCorrCent(X7, 8, 2)  T7["movcorrcent_8_2"]
        m = MovCorr(X7, 10, 1)
        @test size(m) == (size(X7, 1), size(X7, 2) - 1)
        @test all(-1 - 1e-12 .<= m[10:end, :] .<= 1 + 1e-12)
    end

    @testset "argument guards" begin
        @test T7["err_window_zero"] == 1
        @test T7["err_window_big"] == 1
        @test_throws ArgumentError MovAvg(X7, 0)
        @test_throws ArgumentError MovAvg(X7, 1000)
    end
end

@testset "Stats -- censor" begin
    @golden censor(X7, 10) T7["censor_10"]
    @golden censor(X7, 25) T7["censor_25"]

    @testset "censors roughly the requested share" begin
        c = censor(X7, 25)
        for j in 1:size(X7, 2)
            frac = count(isnan, c[:, j]) / size(X7, 1)
            @test 0.4 < frac < 0.6      # 25% from each tail
        end
    end
end

@testset "Stats -- CorrTable" begin
    @golden CorrTable(X7) T7["corrtable"]
    @golden CorrTable(XN7; pairwise = true) T7["corrtable_pair"]

    C = CorrTable(X7)
    @test all(diag(C) .== 1.0)
    @test all(isnan, [C[i, j] for i in 1:size(C, 1) for j in (i+1):size(C, 2)])

    @testset "pairwise uses more data than listwise" begin
        # With NaNs scattered across columns, dropping whole rows loses
        # information that pairwise deletion keeps, so the two differ.
        Cl = CorrTable(XN7)
        Cp = CorrTable(XN7; pairwise = true)
        @test !approxeq(Cl[2, 1], Cp[2, 1]; rtol = 1e-8)
    end
end

@testset "Stats -- PairCorr" begin
    @golden PairCorr(X7) T7["paircorr"]
    @golden PairCorr(X7; absolute = true) T7["paircorr_abs"]

    PC = PairCorr(X7)
    @test size(PC) == (size(X7, 2), 2)
    # Absolute-value averages dominate signed ones.
    PCa = PairCorr(X7; absolute = true)
    @test all(PCa .>= PC .- 1e-12)
end

@testset "Stats -- CrossCorrelation" begin
    cc, lags = CrossCorrelation(X7[:, 1], X7[:, 2:end], 6)
    @golden cc T7["crosscorr"]
    @golden lags T7["crosscorr_lags"]
    @test lags == collect(-6:6)
    @test size(cc) == (13, 3)

    # The centre row is the contemporaneous correlation.
    for ii in 1:3
        @test cc[7, ii] ≈ cor(X7[:, 1], X7[:, ii+1]) rtol = RTOL
    end

    @test_throws DimensionMismatch CrossCorrelation(X7[:, 1], X7[1:10, 2:end], 3)
    @test_throws ArgumentError CrossCorrelation(X7[1:8, 1], X7[1:8, 2:end], 7)
end
