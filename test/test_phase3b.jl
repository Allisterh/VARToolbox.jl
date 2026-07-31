# =======================================================================
# Closes the gap recorded at the end of the Phase 3 log: the
# whereToStart > 1 branch of SignRestrictions.
#
# That branch combines an instrument-identified first column of B with sign
# restrictions on the remaining shocks. It is the only consumer of the
# randn(dy,1) draw site and the only path that reads VAR.sigma_b instead of
# VAR.sigma, so the square-SIGN fixtures in Phase 3 covered none of it.
# =======================================================================

W = fixture("phase3b")

@testset "Phase 3b -- IV + sign restrictions combined" begin
    ENDO = W["ENDO"]
    IVv = W["IV"]
    SIGN = W["SIGN"]

    Vg, og = VARmodel(ENDO, 12, 1)
    Vg = @set Vg.IV = IVv
    og = @set og.nsteps = 12

    IVTOL = 1e-6      # conditioning, as in Phase 3 -- see test_phase3.jl

    @testset "stage 1: instrument populates Biv and sigma_b" begin
        oiv = @set og.ident = "iv"
        _, Vg2 = VARir(Vg, oiv)
        @test approxeq(Vg2.Biv, W["Biv"]; rtol = IVTOL)
        @test approxeq(Vg2.sigma_b, W["sigma_b"]; rtol = IVTOL)
        @test size(Vg2.Biv, 2) == 1
    end

    @testset "stage 2: sign restrictions on the remaining shocks" begin
        oiv = @set og.ident = "iv"
        _, Vg2 = VARir(Vg, oiv)

        @test W["sr_found"] == 1
        @test Int(W["whereToStart"]) == 1 + size(SIGN, 1) - size(SIGN, 2)
        @test Int(W["whereToStart"]) == 2

        osr = @set og.ident = "sign"
        osr = @set osr.sr_hor = 1
        src = DrawSource(normals = [Matrix{Float64}(a) for a in vec(W["sr_draws"])])
        B = SignRestrictions(SIGN, Vg2, osr; src = src)

        @test B !== nothing
        @test approxeq(B, W["sr_B"]; rtol = IVTOL)
    end

    @testset "the randn(dy,1) draw site is actually exercised" begin
        # The whole point of this fixture. In the square-SIGN path every
        # recorded draw is the square rotation matrix; here the orthogonal
        # completion also draws (dy, 1) vectors. If those are absent, this
        # branch was not taken and the test proves nothing.
        sizes = W["sr_draw_sizes"]
        nvar = size(ENDO, 2)
        vec_draws = count(k -> sizes[k, 1] == nvar && sizes[k, 2] == 1, 1:size(sizes, 1))
        sq_draws = count(k -> sizes[k, 1] == sizes[k, 2] > 1, 1:size(sizes, 1))
        @test vec_draws > 0
        @test sq_draws > 0
    end

    @testset "B keeps the instrument-identified first column" begin
        oiv = @set og.ident = "iv"
        _, Vg2 = VARir(Vg, oiv)
        osr = @set og.ident = "sign"
        src = DrawSource(normals = [Matrix{Float64}(a) for a in vec(W["sr_draws"])])
        B = SignRestrictions(SIGN, Vg2, osr; src = src)

        # The first column of B must still be the instrument's, up to the
        # rotation's sign convention -- that is what the two-stage sequence
        # is for.
        @test approxeq(abs.(B[:, 1]), abs.(vec(Vg2.Biv)); rtol = 1e-6)

        # The remaining columns satisfy the imposed signs.
        for j in 1:size(SIGN, 2)
            for i in 1:size(SIGN, 1)
                SIGN[i, j] == 0 && continue
                @test sign(B[i, j+1]) == sign(SIGN[i, j])
            end
        end
    end

    @testset "dimension guards" begin
        oiv = @set og.ident = "iv"
        _, Vg2 = VARir(Vg, oiv)
        osr = @set og.ident = "sign"

        # SIGN wider than tall is rejected outright.
        @test W["err_wide_SIGN"] == 1
        @test_throws DimensionMismatch SignRestrictions(ones(3, 4), Vg2, osr)

        # Biv occupies exactly one column, so SIGN must be exactly one
        # column short of square. A 4x2 SIGN implies two already-identified
        # columns and must be rejected.
        @test W["err_incoherent"] == 1
        @test_throws DimensionMismatch SignRestrictions(
            [0.0 0; 1 0; 1 0; 0 1], Vg2, osr)
    end
end
