# =======================================================================
# Phase 2 golden-master tests: VARoption, VARmodel, ARDLmodel.
#
# Deterministic throughout -- fixtures come from gen_phase2.m, whose input
# data is built from an explicit LCG rather than MATLAB's RNG. Exact
# golden-master at rtol = 1e-10.
# =======================================================================

G = fixture("phase2")
Asyn = G["A"]
Esyn = G["E"]

@testset "Phase 2 -- VARoption defaults" begin
    o = VARoption()
    @test o.nsteps == Int(G["opt_nsteps"])
    @test o.ndraws == Int(G["opt_ndraws"])
    @test o.pctg == G["opt_pctg"]
    @test o.sr_hor == Int(G["opt_sr_hor"])
    @test o.sr_rot == Int(G["opt_sr_rot"])
    @test o.sr_draw == Int(G["opt_sr_draw"])
    @test o.sr_mod == Int(G["opt_sr_mod"])
    @test o.mult == Int(G["opt_mult"])
    @test o.quality == Int(G["opt_quality"])
    @test o.impact == Int(G["opt_impact"])
    @test o.shut == Int(G["opt_shut"])
    @test o.pick == Int(G["opt_pick"])
    @test o.datestype == Int(G["opt_datestype"])
    @test o.suptitle == Int(G["opt_suptitle"])
    @test o.ident == G["opt_ident"]
    @test o.recurs == G["opt_recurs"]
    @test o.method == G["opt_method"]
    @test o.frequency == G["opt_frequency"]
    @test o.FigSize == flat(G["opt_FigSize"])

    @testset "immutability (D3)" begin
        # VAROptions must be immutable: SR.m and SignRestrictions.m mutate
        # their local copy and rely on MATLAB's by-value semantics to keep
        # the caller's options clean.
        @test !ismutable(o)
        o2 = @set o.nsteps = 24
        @test o2.nsteps == 24
        @test o.nsteps == 40      # original untouched
    end
end

@testset "Phase 2 -- VARmodel (synthetic)" begin
    cfg = [(2, 1), (1, 1), (3, 1), (2, 0), (2, 2), (2, 3)]
    @test length(cfg) == Int(G["nconfig"])

    for (k, (nlag, cst)) in enumerate(cfg)
        @testset "nlag=$nlag const=$cst" begin
            V, opt = VARmodel(Asyn, nlag, cst)
            p = "v$(k)_"
            @test V.nlag == Int(G[p*"nlag"])
            @test V.constant == Int(G[p*"const"])
            @test V.nobs == Int(G[p*"nobs"])
            @test V.nvar == Int(G[p*"nvar"])
            @test V.ncoeff == Int(G[p*"ncoeff"])
            @test V.ntotcoeff == Int(G[p*"ntotcoeff"])
            @golden V.Ft    G[p*"Ft"]
            @golden V.F     G[p*"F"]
            @golden V.sigma G[p*"sigma"]
            @golden V.resid G[p*"resid"]
            @golden V.X     G[p*"X"]
            @golden V.Y     G[p*"Y"]
            @golden V.Fcomp G[p*"Fcomp"]
            @test V.maxEig ≈ G[p*"maxEig"] rtol = RTOL

            # Identification is late-bound: nothing filled in yet.
            @test isempty(V.B)
            @test isempty(V.Biv)
            @test isempty(V.PSI)
            @test isempty(V.IV)

            # The 16 eval-built VAR.eqN fields become a typed vector.
            @test length(V.eq) == V.nvar
            for j in 1:V.nvar
                q = "$(p)eq$(j)_"
                @golden V.eq[j].beta  G[q*"beta"]
                @golden V.eq[j].tstat G[q*"tstat"]
                @golden V.eq[j].bstd  G[q*"bstd"]
                # tprob: v3.1 used nobs - ncoeff (lags only) degrees of
                # freedom; the port follows v4.0 (nobs - ntotcoeff) and is
                # checked against the v4.0 fixture in test_phase9.jl.
                @golden V.eq[j].resid G[q*"resid"]
                @golden V.eq[j].yhat  G[q*"yhat"]
                @golden V.eq[j].y     G[q*"y"]
                @test V.eq[j].rsqr ≈ G[q*"rsqr"] rtol = RTOL
                @test V.eq[j].rbar ≈ G[q*"rbar"] rtol = RTOL
                @test V.eq[j].sige ≈ G[q*"sige"] rtol = RTOL
                @test V.eq[j].dw   ≈ G[q*"dw"]   rtol = RTOL
            end
        end
    end

    @testset "const_ defaults to 1, not 0" begin
        # VARmodel.m's docstring says [dflt = 0]; the code uses 1. The code
        # is authoritative -- this is one of the two agreed doc-fix carve-outs.
        V1, _ = VARmodel(Asyn, 2)
        V2, _ = VARmodel(Asyn, 2, 1)
        @test V1.constant == 1
        @test V1.ntotcoeff == V2.ntotcoeff
        @golden V1.Ft V2.Ft
    end

    @testset "second return value is default options" begin
        _, opt = VARmodel(Asyn, 2)
        @test opt == VARoption()
    end

    @test_throws ArgumentError VARmodel(Asyn, 0)
end

@testset "Phase 2 -- VARmodel with EXOG" begin
    # All three nlag vs nlag_ex branches: equal, longer, shorter.
    ex = [(2, 2), (3, 1), (1, 3)]
    @test length(ex) == Int(G["nexog"])

    for (k, (nlag, nlx)) in enumerate(ex)
        @testset "nlag=$nlag nlag_ex=$nlx" begin
            V, _ = VARmodel(Asyn, nlag, 1, Esyn, nlx)
            p = "x$(k)_"
            @test V.nlag == Int(G[p*"nlag"])
            @test V.nlag_ex == Int(G[p*"nlag_ex"])
            @test V.nobs == Int(G[p*"nobs"])
            @test V.nvar_ex == Int(G[p*"nvar_ex"])
            @test V.ncoeff == Int(G[p*"ncoeff"])
            @test V.ntotcoeff == Int(G[p*"ntotcoeff"])
            @golden V.Ft    G[p*"Ft"]
            @golden V.F     G[p*"F"]
            @golden V.sigma G[p*"sigma"]
            @golden V.resid G[p*"resid"]
            @golden V.X     G[p*"X"]
            @golden V.Y     G[p*"Y"]
            @golden V.X_EX  G[p*"X_EX"]
            @golden V.Fcomp G[p*"Fcomp"]
            @test V.maxEig ≈ G[p*"maxEig"] rtol = RTOL
        end
    end

    @test_throws DimensionMismatch VARmodel(Asyn, 2, 1, Esyn[1:10, :], 1)
end

@testset "Phase 2 -- VARmodel on real data (SW2001)" begin
    # Exactly the configuration GO_SW2001.m uses. Validating here means a
    # Phase 6 end-to-end discrepancy can be localised to identification
    # rather than estimation.
    SW = G["SW"]
    V, _ = VARmodel(SW, 4, 1)
    @test V.nobs == Int(G["sw_nobs"])
    @test V.nvar == Int(G["sw_nvar"])
    @test V.ntotcoeff == Int(G["sw_ntotcoeff"])
    @golden V.Ft    G["sw_Ft"]
    @golden V.F     G["sw_F"]
    @golden V.sigma G["sw_sigma"]
    @golden V.resid G["sw_resid"]
    @golden V.Fcomp G["sw_Fcomp"]
    @golden V.X     G["sw_X"]
    @golden V.Y     G["sw_Y"]
    @test V.maxEig ≈ G["sw_maxEig"] rtol = RTOL

    # A VAR fitted to real macro data should be stable.
    @test V.maxEig < 1.0

    for j in 1:V.nvar
        q = "sw_eq$(j)_"
        @golden V.eq[j].beta  G[q*"beta"]
        @golden V.eq[j].tstat G[q*"tstat"]
        # tprob checked against v4.0 in test_phase9.jl (see above).
        @test V.eq[j].rsqr ≈ G[q*"rsqr"] rtol = RTOL
        @test V.eq[j].dw   ≈ G[q*"dw"]   rtol = RTOL
        @test V.eq[j].sige ≈ G[q*"sige"] rtol = RTOL
    end
end

@testset "Phase 2 -- ARDLmodel" begin
    y1 = Asyn[:, 1:1]

    @testset "no EXOG" begin
        D = ARDLmodel(y1, 2, 1)
        @test D.nobs == Int(G["ardl_nobs"])
        @test D.nvar == Int(G["ardl_nvar"])
        @test D.ncoeff == Int(G["ardl_ncoeff"])
        @golden D.beta    G["ardl_beta"]
        @golden D.yhat    G["ardl_yhat"]
        @golden D.resid   G["ardl_resid"]
        @golden D.bstd    G["ardl_bstd"]
        # bint / tprob: v3.1 used nobse and nobs degrees of freedom; the
        # port follows v4.0 (nobse - nvar), checked in test_phase9.jl.
        @golden D.tstat   G["ardl_tstat"]
        @golden D.sigbeta G["ardl_sigbeta"]
        @golden D.x       G["ardl_x"]
        @golden D.y       G["ardl_y"]
        @test D.sige ≈ G["ardl_sige"] rtol = RTOL
        @test D.rsqr ≈ G["ardl_rsqr"] rtol = RTOL
        @test D.rbar ≈ G["ardl_rbar"] rtol = RTOL
        @test D.dw   ≈ G["ardl_dw"]   rtol = RTOL
        @test D.F    ≈ G["ardl_F"]    rtol = RTOL
    end

    @testset "with EXOG" begin
        D = ARDLmodel(y1, 2, 1, Esyn, 1)
        @test D.nobs == Int(G["ardl2_nobs"])
        @test D.nvar == Int(G["ardl2_nvar"])
        @golden D.beta  G["ardl2_beta"]
        @golden D.resid G["ardl2_resid"]
        @golden D.tstat G["ardl2_tstat"]
        @golden D.x     G["ardl2_x"]
        @golden D.y     G["ardl2_y"]
        @test D.rsqr ≈ G["ardl2_rsqr"] rtol = RTOL
        @test D.dw   ≈ G["ardl2_dw"]   rtol = RTOL
    end

    @testset "no constant" begin
        D = ARDLmodel(y1, 2, 0)
        @golden D.beta G["ardl3_beta"]
        @test D.nvar == Int(G["ardl3_nvar"])
        @test D.rsqr ≈ G["ardl3_rsqr"] rtol = RTOL
        @test isnan(D.F)
    end

    @testset "rejects multivariate input" begin
        @test_throws ArgumentError ARDLmodel(Asyn, 2, 1)
    end
end
