# =======================================================================
# Phase 3 golden-master tests: VARir, VARvd, VARhd, SignRestrictions, SR,
# VARdrawpost.
#
# The sign-restriction tests are EXACT golden-master checks, not
# distributional ones: gen_phase3.m records the randn sequence MATLAB
# consumed and this file replays it through a DrawSource (M2).
# =======================================================================

H = fixture("phase3")

"MATLAB cell array of recorded draws -> Vector{Matrix{Float64}}"
drawpool(x) = [Matrix{Float64}(a) for a in vec(x)]

SWdata = H["SW"]
Vsw, osw = VARmodel(SWdata, 4, 1)
osw = @set osw.nsteps = 12

@testset "Phase 3 -- VARir" begin
    @testset "short / long" begin
        for nm in ("short", "long")
            o = @set osw.ident = nm
            IR, Vo = VARir(Vsw, o)
            @golden IR    H["ir_$nm"]
            @golden Vo.B  H["B_$nm"]
            @golden Vo.PSI H["PSI_$nm"]
            @golden Vo.Fp  H["Fp_$nm"]
            @test size(IR) == (12, Vsw.nvar, Vsw.nvar)
        end
    end

    @testset "IR layout is (nsteps, nvar, nshocks)" begin
        o = @set osw.ident = "short"
        IR, _ = VARir(Vsw, o)
        # Impact response to shock 1 is column 1 of B.
        _, Vo = VARir(Vsw, o)
        @test IR[1, :, 1] ≈ Vo.B[:, 1] rtol = RTOL
    end

    @testset "recurs = comp matches wold" begin
        o = @set osw.ident = "short"
        o = @set o.recurs = "comp"
        @golden VARir(Vsw, o)[1] H["ir_short_comp"]
        # The two recursions are algebraically the same object.
        ow = @set osw.ident = "short"
        @test approxeq(VARir(Vsw, o)[1], VARir(Vsw, ow)[1]; rtol = 1e-8)
    end

    @testset "impact = 1 (unit shock)" begin
        o = @set osw.ident = "short"
        o = @set o.impact = 1
        IR, Vo = VARir(Vsw, o)
        @golden IR H["ir_short_unit"]
        # A unit shock moves variable mm by exactly 1 on impact.
        for mm in 1:Vsw.nvar
            @test IR[1, mm, mm] ≈ 1.0 rtol = 1e-10
        end
    end

    @testset "shut" begin
        o = @set osw.ident = "short"
        o = @set o.shut = 2
        IR, _ = VARir(Vsw, o)
        @golden IR H["ir_short_shut2"]
        @test all(IR[1, 2, :] .== 0)
    end

    @testset "ident = sign without B raises" begin
        @test H["ir_sign_no_B_errors"] == 1
        o = @set osw.ident = "sign"
        @test_throws ErrorException VARir(Vsw, o)
    end

    @testset "unknown ident raises" begin
        o = @set osw.ident = "nonsense"
        @test_throws ErrorException VARir(Vsw, o)
    end
end

@testset "Phase 3 -- VARvd" begin
    for nm in ("short", "long")
        o = @set osw.ident = nm
        VD, _ = VARvd(Vsw, o)
        @golden VD H["vd_$nm"]
    end

    @testset "VD layout is (nsteps, nshocks, nvar) -- NOT IR's layout" begin
        o = @set osw.ident = "short"
        VD, _ = VARvd(Vsw, o)
        # Contributions sum to 100 across the SECOND index.
        @test sum(VD[1, :, 1]) ≈ H["vd_sum_dim2"] rtol = RTOL
        @test sum(VD[1, :, 1]) ≈ 100.0 rtol = 1e-10
        @test sum(VD[1, 1, :]) ≈ H["vd_sum_dim3"] rtol = RTOL
        @test !isapprox(sum(VD[1, 1, :]), 100.0; rtol = 1e-6)
        for h in 1:12, k in 1:Vsw.nvar
            @test sum(VD[h, :, k]) ≈ 100.0 rtol = 1e-9
        end
    end

    @testset "rejects iv and exog" begin
        @test H["vd_rejects_iv"] == 1
        @test H["vd_rejects_exog"] == 1
        @test_throws ErrorException VARvd(Vsw, @set osw.ident = "iv")
        @test_throws ErrorException VARvd(Vsw, @set osw.ident = "exog")
    end
end

@testset "Phase 3 -- VARhd" begin
    for nm in ("short", "long")
        o = @set osw.ident = nm
        HD, _ = VARhd(Vsw, o)
        @golden HD.shock    H["hd_$(nm)_shock"]
        @golden HD.init     H["hd_$(nm)_init"]
        @golden HD.constant H["hd_$(nm)_const"]
        @golden HD.trend    H["hd_$(nm)_trend"]
        @golden HD.trend2   H["hd_$(nm)_trend2"]
        @golden HD.endo     H["hd_$(nm)_endo"]
    end

    @testset "NaN padding" begin
        o = @set osw.ident = "short"
        HD, _ = VARhd(Vsw, o)
        nlag = Vsw.nlag
        # init pads with nlag-1; everything else with nlag. Upstream quirk.
        @test all(isnan, HD.init[1:(nlag-1), :])
        @test !any(isnan, HD.init[nlag:end, :])
        @test all(isnan, HD.constant[1:nlag, :])
        @test !any(isnan, HD.constant[(nlag+1):end, :])
        @test size(HD.endo, 1) == Vsw.nobs + nlag
    end

    @testset "components reconstruct the data" begin
        # With no exogenous variables, endo = init + const + trend + trend2
        # + sum of shock contributions, and that should reproduce Y.
        o = @set osw.ident = "short"
        HD, _ = VARhd(Vsw, o)
        nlag = Vsw.nlag
        recon = HD.endo[(nlag+1):end, :]
        @test approxeq(recon, Vsw.Y; rtol = 1e-8)
    end

    @testset "rejects iv and exog" begin
        @test H["hd_rejects_iv"] == 1
        @test H["hd_rejects_exog"] == 1
        @test_throws ErrorException VARhd(Vsw, @set osw.ident = "iv")
        @test_throws ErrorException VARhd(Vsw, @set osw.ident = "exog")
    end
end

@testset "Phase 3 -- ident = exog" begin
    EX = reshape(flat(H["EX"]), :, 1)
    Vx, ox = VARmodel(SWdata, 4, 1, EX, 0)
    ox = @set ox.nsteps = 12
    @test Vx.nvar_ex == Int(H["vx_nvar_ex"])
    @test Vx.ntotcoeff == Int(H["vx_ntotcoeff"])

    IRx, Vxo = VARir(Vx, @set ox.ident = "exog")
    @golden IRx    H["ir_exog"]
    @golden Vxo.B  H["B_exog"]
    # exog fills only the first column of B.
    @test all(Vxo.B[:, 2:end] .== 0)

    HDx, _ = VARhd(Vx, @set ox.ident = "short")
    @golden HDx.exo   H["hdx_exo"]
    @golden HDx.endo  H["hdx_endo"]
    @golden HDx.shock H["hdx_shock"]
end

@testset "Phase 3 -- ident = iv (Gertler-Karadi)" begin
    ENDO = H["gk_ENDO"]
    IVv = H["gk_IV"]
    Vg, og = VARmodel(ENDO, 12, 1)
    Vg = @set Vg.IV = IVv
    og = @set og.nsteps = 12
    og = @set og.ident = "iv"

    # ---- tolerance justification (required by the validation policy) ----
    #
    # This block runs at rtol = 1e-6 rather than 1e-10. The reason is
    # conditioning, not a porting discrepancy.
    #
    # GO_GK2015.m fits 12 lags of 4 variables plus a constant: 49
    # coefficients on 384 observations, from persistent macro series. The
    # resulting design has cond(X) ~ 4.7e4, so the normal-equations matrix
    # X'X has cond ~ 2.2e9. Solving through X'X therefore admits a relative
    # error of order cond(X'X)*eps ~ 5e-7 regardless of implementation.
    #
    # Measured difference against MATLAB is 6.5e-8 in norm-relative terms
    # (3.1e-6 absolute on coefficients up to 32.1) -- comfortably inside
    # what the algorithm supports, and an order of magnitude better than the
    # bound. Both answers are equally valid solutions of the same
    # ill-conditioned system.
    #
    # The assertion below makes this justification self-checking: if the
    # design ever becomes well conditioned, the loosened tolerance is no
    # longer warranted and this test should be tightened.
    IVTOL = 1e-6
    @test cond(Vg.X' * Vg.X) > 1e8
    @test cond(Vg.X' * Vg.X) * eps() > 1e-8   # 1e-10 is genuinely unreachable here

    @test approxeq(Vg.sigma, H["gk_sigma"]; rtol = IVTOL)
    @test approxeq(Vg.Ft, H["gk_Ft"]; rtol = IVTOL)

    IRiv, Vgo = VARir(Vg, og)
    @test approxeq(IRiv, H["ir_iv"]; rtol = IVTOL)
    @test approxeq(Vgo.B, H["B_iv"]; rtol = IVTOL)
    @test approxeq(Vgo.Biv, H["Biv"]; rtol = IVTOL)
    @test approxeq(Vgo.sigma_b, H["sigma_b"]; rtol = IVTOL)
    @test approxeq(Vgo.FirstStage.beta, H["gk_fs_beta"]; rtol = IVTOL)
    @test approxeq(Vgo.FirstStage.yhat, H["gk_fs_yhat"]; rtol = IVTOL)
    @test Vgo.FirstStage.rsqr ≈ H["gk_fs_rsqr"] rtol = IVTOL

    # iv fills only the first column of B.
    @test all(Vgo.B[:, 2:end] .== 0)

    @testset "variable order is load-bearing" begin
        # The instrumented variable must be ordered FIRST. Swapping it out
        # of position 1 must change the answer -- proving the convention is
        # real rather than incidental.
        ENDO2 = ENDO[:, [2, 1, 3, 4]]
        Vg2, og2 = VARmodel(ENDO2, 12, 1)
        Vg2 = @set Vg2.IV = IVv
        og2 = @set og2.nsteps = 12
        og2 = @set og2.ident = "iv"
        IRsw, _ = VARir(Vg2, og2)
        @test approxeq(IRsw, H["ir_iv_swapped"]; rtol = IVTOL)
        @test !approxeq(IRsw, IRiv; rtol = 1e-3)
    end

    @testset "missing IV raises" begin
        Vno, ono = VARmodel(ENDO, 12, 1)
        ono = @set ono.ident = "iv"
        @test_throws ErrorException VARir(Vno, ono)
    end
end

@testset "Phase 3 -- SignRestrictions (M2 exact replay)" begin
    SIGN = H["SIGN"]

    @testset "sr_hor = 1" begin
        src = DrawSource(normals = drawpool(H["sr_draws"]))
        o = @set osw.ident = "sign"
        o = @set o.sr_hor = 1
        B = SignRestrictions(SIGN, Vsw, o; src = src)
        @test B !== nothing
        @golden B H["sr_B"]
        # Every recorded draw was consumed -- confirms the two
        # implementations took the same number of rotations.
        @test src.i == Int(H["sr_ndraws_used"])
    end

    @testset "sr_hor = 4 (checks IRFs, not just impact)" begin
        src = DrawSource(normals = drawpool(H["sr2_draws"]))
        o = @set osw.ident = "sign"
        o = @set o.sr_hor = 4
        B = SignRestrictions(SIGN, Vsw, o; src = src)
        @test B !== nothing
        @golden B H["sr2_B"]
    end

    @testset "B satisfies the sign pattern" begin
        src = DrawSource(normals = drawpool(H["sr_draws"]))
        o = @set osw.ident = "sign"
        B = SignRestrictions(SIGN, Vsw, o; src = src)
        for i in 1:size(SIGN, 1)
            SIGN[i, 1] != 0 && @test sign(B[i, 1]) == sign(SIGN[i, 1])
        end
        # B is a valid factorisation of sigma.
        @test approxeq(B * B', Vsw.sigma; rtol = 1e-9)
    end

    @testset "exhausted rotation budget returns nothing, silently" begin
        # Upstream's error path is commented out -- it sets flag=1 and
        # returns []. Impossible restrictions must therefore yield `nothing`
        # and NOT raise.
        o = @set osw.ident = "sign"
        o = @set o.sr_rot = 5
        impossible = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
        B = SignRestrictions(impossible, Vsw, o; src = DrawSource(rng = MersenneTwister(1)))
        @test B === nothing
    end
end

@testset "Phase 3 -- SR end to end (sr_mod = 0, exact)" begin
    SIGN = H["SIGN"]
    src = DrawSource(normals = drawpool(H["SR_draws"]))
    o = @set osw.ident = "sign"
    o = @set o.ndraws = Int(H["SR_ndraws"])
    o = @set o.sr_mod = 0
    o = @set o.nsteps = Int(H["SR_nsteps"])

    out = SR(Vsw, SIGN, o; src = src)

    @golden out.Ball  H["SR_Ball"]
    @golden out.Bmed  H["SR_Bmed"]
    @golden out.B     H["SR_B"]
    @golden out.IRall H["SR_IRall"]
    @golden out.VDall H["SR_VDall"]
    @golden out.IRmed H["SR_IRmed"]
    @golden out.IRinf H["SR_IRinf"]
    @golden out.IRsup H["SR_IRsup"]
    @golden out.VDmed H["SR_VDmed"]
    @golden out.VDinf H["SR_VDinf"]
    @golden out.VDsup H["SR_VDsup"]
    @golden out.IR    H["SR_IR"]
    @golden out.VD    H["SR_VD"]
    @golden out.HD.shock H["SR_HD_shock"]
    @golden out.HD.endo  H["SR_HD_endo"]
    @golden out.HDall[:endo] H["SR_HDall_endo"]

    @testset "upstream off-by-one between B and IR" begin
        # SR.m labels each stored draw by the accepted count BEFORE
        # incrementing it, so the model it reads back at the end is the one
        # AFTER the draw that B came from. Verified in MATLAB: B was
        # accepted draw 24 while IR matched draw 25 exactly.
        selB = Int(H["SR_selB"])
        selIR = Int(H["SR_selIR"])
        @test H["SR_selIRresid"] < 1e-12       # IR really is some draw, exactly
        @test selIR == selB + 1                # and it is the NEXT one
        @test selB != selIR

        # The port reproduces the mismatch by default...
        @test approxeq(out.B, H["SR_B"])
        @test !approxeq(out.IR, VARir((@set Vsw.B = out.B), (@set osw.ident = "sign"))[1];
                        rtol = 1e-6)

        # ...and fix_draw_index=true makes IR correspond to B.
        src2 = DrawSource(normals = drawpool(H["SR_draws"]))
        fixed = SR(Vsw, SIGN, o; src = src2, fix_draw_index = true)
        @test approxeq(fixed.B, out.B)         # same selected B
        expected, _ = VARir((@set Vsw.B = fixed.B), (@set osw.ident = "sign"))
        @test approxeq(fixed.IR, expected[1:Int(H["SR_nsteps"]), :, :]; rtol = 1e-9)
    end

    @testset "bands bracket the median" begin
        @test all(out.IRinf .<= out.IRmed .+ 1e-12)
        @test all(out.IRmed .<= out.IRsup .+ 1e-12)
    end

    @testset "exhausted draw budget raises" begin
        o2 = @set o.sr_draw = 3
        o2 = @set o2.ndraws = 1000
        impossible = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
        @test_throws ErrorException SR(Vsw, impossible, o2;
                                       rng = MersenneTwister(2))
    end
end

@testset "Phase 3 -- VARdrawpost" begin
    # The draw itself uses wishrnd/mvnrnd and cannot match value-for-value.
    # Inject MATLAB's draw and validate everything downstream of it exactly.
    sd = H["vdp_sigma_draw"]
    fd = H["vdp_Ft_draw"]
    _, Ft_draw, F_draw, Fcomp_draw =
        VARdrawpost(Vsw; sigma_draw = sd, Ft_draw = fd)
    @golden Ft_draw    H["vdp_Ft_draw"]
    @golden F_draw     H["vdp_F_draw"]
    @golden Fcomp_draw H["vdp_Fcomp_draw"]

    @testset "sampling path runs and is well formed" begin
        rng = MersenneTwister(3)
        sdr, ftr, fdr, fcr = VARdrawpost(Vsw; rng = rng)
        @test size(sdr) == size(Vsw.sigma)
        @test issymmetric(Symmetric(sdr))
        @test all(diag(sdr) .> 0)
        @test size(ftr) == size(Vsw.Ft)
        @test size(fcr) == size(Vsw.Fcomp)
    end
end
