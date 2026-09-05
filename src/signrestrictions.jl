# =======================================================================
# signrestrictions.jl
#
# Port of v3dot1/VAR/SignRestrictions.m, SR.m and VARdrawpost.m.
#
# This is where M2 (injectable randomness) lives. Every stochastic step
# takes its draws through a `DrawSource`, so MATLAB's own draws can be fed
# in and the whole sign-restriction pipeline validated as an exact
# golden-master rather than a distributional comparison.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

"""
    DrawSource(; rng=Random.default_rng(), normals=nothing)

Source of random draws for the stochastic routines (M2).

With `normals = nothing` it samples from `rng`, which is the normal
production path. With `normals` set to a vector of pre-generated arrays it
returns them in order instead — that is how MATLAB's draws are injected so
that bootstrap and rotation results can be compared exactly rather than
distributionally.

Exhausting an injected pool is an error, not a silent fallback to sampling:
a test that quietly stopped being a golden-master check would be worse than
no test.
"""
mutable struct DrawSource
    rng::AbstractRNG
    normals::Union{Nothing,Vector{<:AbstractArray}}
    i::Int
    uniforms::Union{Nothing,Vector{<:AbstractArray}}
    j::Int
end
DrawSource(; rng::AbstractRNG = Random.default_rng(), normals = nothing,
           uniforms = nothing) = DrawSource(rng, normals, 0, uniforms, 0)

function next_normal!(s::DrawSource, dims::Integer...)
    if s.normals === nothing
        return randn(s.rng, dims...)
    end
    s.i += 1
    s.i > length(s.normals) &&
        error("DrawSource: injected pool exhausted after $(s.i - 1) draws; " *
              "generate a longer pool in the fixture generator")
    a = s.normals[s.i]
    size(a) == dims ||
        error("DrawSource: injected draw $(s.i) has size $(size(a)), expected $(dims)")
    return copy(a)
end

# -----------------------------------------------------------------------
# VARdrawpost
# -----------------------------------------------------------------------

"""
    VARdrawpost(VAR; rng=Random.default_rng(), sigma_draw=nothing, Ft_draw=nothing)
        -> (sigma_draw, Ft_draw, F_draw, Fcomp_draw)

One draw of `(sigma, F)` from the posterior, used by [`SR`](@ref) when
`VARopt.sr_mod == 1` to add model uncertainty on top of identification
uncertainty.

Maps MATLAB's `wishrnd`/`mvnrnd` to `Distributions.jl`. Note the
parameterisation: MATLAB's `wishrnd(S, df)` takes the **scale** matrix, and
`Wishart(df, S)` in `Distributions.jl` matches it (both have mean `df*S`).

Pass `sigma_draw` and/or `Ft_draw` explicitly to bypass sampling — used to
inject MATLAB's draws for golden-master testing of everything downstream of
the draw itself.

!!! note "This is the one place the Statistics Toolbox is genuinely required"
    `wishrnd` and `mvnrnd` are not in base MATLAB. Fixtures for this path
    were generated on an installation with the Statistics and Machine
    Learning Toolbox licensed (see PROVENANCE.md).
"""
function VARdrawpost(VAR::VARModel; rng::AbstractRNG = Random.default_rng(),
                     sigma_draw = nothing, Ft_draw = nothing)
    nobs, nvar, nlag = VAR.nobs, VAR.nvar, VAR.nlag
    X = VAR.X
    constant = VAR.constant

    if sigma_draw === nothing
        inv_sigma_hat = inv(VAR.sigma)
        # Symmetrise before handing to Wishart: inv() of a symmetric matrix
        # is only symmetric up to rounding, and Wishart requires an exactly
        # symmetric (positive definite) scale matrix.
        S = Matrix(Symmetric(inv_sigma_hat ./ nobs))
        inv_sigma_draw = rand(rng, Wishart(float(nobs), S))
        sigma_draw = inv(Matrix(inv_sigma_draw))
    end

    if Ft_draw === nothing
        aux1 = inv(transpose(X) * X)
        aux2 = kron(sigma_draw, aux1)
        aux2 = (aux2 + transpose(aux2)) ./ 2      # force symmetry (upstream)
        Ftdraw = rand(rng, MvNormal(vec(VAR.Ft), Symmetric(aux2)))
        Ft_draw = reshape(Ftdraw, :, nvar)
    end

    F_draw = Matrix(transpose(Ft_draw))
    k = nvar * (nlag - 1)
    top = F_draw[:, (1+constant):(nvar*nlag+constant)]
    Fcomp_draw = k == 0 ? top :
                 vcat(top, hcat(Matrix{Float64}(I, k, k), zeros(k, nvar)))
    return sigma_draw, Ft_draw, F_draw, Fcomp_draw
end

# -----------------------------------------------------------------------
# SignRestrictions
# -----------------------------------------------------------------------

"""
    SignRestrictions(SIGN, VAR, VARopt; src=DrawSource()) -> B or nothing

Draw rotations of the Cholesky factor until one satisfies the sign pattern
in `SIGN`, or until `VARopt.sr_rot` rotations have been tried.

`SIGN` is `nvar × nshocks`; entry `(i, j)` is `+1`, `-1` or `0` for the
required sign of variable `i`'s response to shock `j`. When `SIGN` has fewer
columns than rows, the leading shocks are taken as already identified (see
[`SR`](@ref) for combining with an external instrument).

Returns `nothing` when the rotation budget is exhausted. Exactly
`VARopt.sr_rot` rotations are attempted.

!!! note "Rotation budget off-by-one (upstream v4.0 fix)"
    VAR Toolbox 3.1 incremented its rotation counter *before* testing it
    against `sr_rot`, so at most `sr_rot - 1` rotations were tried and
    `sr_rot = 1` returned `[]` unconditionally — `SR` then treated every
    call as a rejection and cycled until the `sr_draw` budget aborted with
    a misleading "restrictions too tight" message. v4.0 tests the budget
    first (`ChangeLog.md`, 2026-08-16); this port follows v4.0. The draw
    sequence is unchanged whenever the restrictions are met within the
    budget, so recorded MATLAB draws still replay exactly.

!!! warning "Budget exhaustion is silent upstream"
    `SignRestrictions.m` sets `flag = 1` and returns `[]` when it exceeds
    `VARopt.sr_rot`; its `error(...)` call is **commented out** in the
    source. So a caller that does not check for an empty return gets no
    warning at all. This port returns `nothing` and does **not** raise,
    matching that behaviour — [`SR`](@ref) is the layer that turns repeated
    failures into a hard error, via its separate `sr_draw` budget.

    Note this is a *different* budget from `sr_draw`: `sr_rot` bounds
    rotations within one call, `sr_draw` bounds total calls across `SR`.

If `VARopt.sr_hor > 1` the signs are checked on the impulse responses over
`sr_hor` horizons rather than on impact. Note that this path calls
[`VARir`](@ref) with `nsteps` temporarily set to `sr_hor` on a local copy of
the options — one of the mutations that makes [`VAROptions`](@ref)'s
immutability load-bearing.
"""
function SignRestrictions(SIGN::AbstractMatrix, VAR::VARModel, VARopt::VAROptions;
                          src::DrawSource = DrawSource())
    dy, ds = size(SIGN)

    if isempty(VAR.Biv)
        if dy == ds
            sigma = VAR.sigma
            Biv = zeros(dy, 0)
        elseif dy > ds
            sigma = VAR.sigma
            aux = Matrix(cholesky(Symmetric(VAR.sigma)).U)'
            Biv = aux[:, (dy-ds):(dy-ds)]
        else
            throw(DimensionMismatch("Matrix SIGN has the wrong dimension"))
        end
    else
        sigma = VAR.sigma_b
        Biv = VAR.Biv
    end

    whereToStart = 1 + dy - ds
    size(Biv, 2) == whereToStart - 1 ||
        throw(DimensionMismatch("Biv and SIGN must have coherent sizes"))

    orderIndices = collect(1:dy)
    counter = 0
    flag = 0
    termaa = zeros(dy, dy)

    while true
        # v4.0: test the budget BEFORE building the rotation, so exactly
        # sr_rot rotations are attempted (v3.1 attempted sr_rot - 1).
        if counter >= VARopt.sr_rot
            flag = 1
            break
        end
        counter += 1

        if whereToStart > 1
            C = Matrix(cholesky(Symmetric(sigma)).L)
            q = C \ Biv
            for _ in whereToStart:dy
                r = next_normal!(src, dy, 1)
                w = (Matrix{Float64}(I, dy, dy) - q * transpose(q)) * r
                q = hcat(q, w ./ norm(w))
            end
            startingMat = C * q
        else
            startingMat = Matrix(cholesky(Symmetric(sigma)).L)
        end

        termaa = copy(startingMat)
        TermA = 0
        rotMat = Matrix{Float64}(I, dy, dy)
        nrot = length(whereToStart:dy)
        rotMat[whereToStart:end, whereToStart:end] = getqr(next_normal!(src, nrot, nrot))
        terma = termaa * rotMat
        termaa = copy(terma)

        aux_irf = nothing
        if VARopt.sr_hor > 1
            # Local copy only -- the caller's options must not see this.
            opt_h = @set VARopt.nsteps = VARopt.sr_hor
            VARh = @set VAR.B = termaa
            aux_irf, _ = VARir(VARh, opt_h)
        end

        for ii in 1:ds
            for jj in whereToStart:dy
                isfinite(terma[1, jj]) || continue
                CHECK = VARopt.sr_hor > 1 ? transpose(aux_irf[:, :, jj]) : terma[:, jj:jj]
                if !any(CHECK .* SIGN[:, ii] .< 0)
                    TermA += 1
                    orderIndices[whereToStart-1+ii] = jj
                    terma[:, jj] .= NaN
                    break
                elseif !any(.-CHECK .* SIGN[:, ii] .< 0)
                    TermA += 1
                    terma[:, jj] .= NaN
                    termaa[:, jj] .= .-termaa[:, jj]
                    orderIndices[ii+whereToStart-1] = jj
                    break
                end
            end
        end

        TermA == ds && break
    end

    return flag == 0 ? termaa[:, orderIndices] : nothing
end

# -----------------------------------------------------------------------
# SR
# -----------------------------------------------------------------------

"""
    SROutput

Result of [`SR`](@ref).

`IRall`/`VDall` are `(nsteps × · × · × ndraws)`; `Ball` is
`(nvar × nvar × ndraws)`. The `med`/`inf`/`sup` fields are the pointwise
median and `VARopt.pctg` band across draws.

`IR`, `VD` and `HD` are computed from the same accepted draw as `B` unless
`SR` was called with `fix_draw_index = false` (the VAR Toolbox 3.1
behaviour); see [`SR`](@ref).
"""
Base.@kwdef struct SROutput
    IRall::Array{Float64,4}
    VDall::Array{Float64,4}
    Ball::Array{Float64,3}
    HDall::Dict{Symbol,Array{Float64}}
    Bmed::Matrix{Float64}
    B::Matrix{Float64}
    IRmed::Array{Float64,3}
    IRinf::Array{Float64,3}
    IRsup::Array{Float64,3}
    VDmed::Array{Float64,3}
    VDinf::Array{Float64,3}
    VDsup::Array{Float64,3}
    IR::Array{Float64,3}
    VD::Array{Float64,3}
    HD::HDOutput
end

"""
    SR(VAR, SIGN, VARopt; rng=Random.default_rng(), src=nothing,
       post_draws=nothing, fix_draw_index=true) -> SROutput

Sign-restriction identification. Collects `VARopt.ndraws` accepted rotations
and returns impulse responses, variance decompositions and historical
decompositions across them.

This **bypasses [`VARir`](@ref)'s `ident` dispatch** — call `SR` rather than
setting `ident = "sign"` by hand.

When `VARopt.sr_mod == 1` (the default) each draw also redraws `F` and
`sigma` from the posterior via [`VARdrawpost`](@ref), so the bands reflect
model *and* identification uncertainty. Set `sr_mod = 0` to hold the
estimated model fixed.

!!! note "Posterior draws carry their own residuals (upstream v4.0 fix)"
    Under `sr_mod == 1` the residuals of each drawn model are recomputed as
    `Y - X*Ft_draw`, so the historical decomposition of a draw uses the
    innovations implied by that draw's coefficients and `HD.endo`
    reconstructs the data. VAR Toolbox 3.1 replaced `Ft`, `F`, `Fcomp` and
    `sigma` but left `resid` at the OLS point estimate, mixing drawn
    dynamics with point-estimate innovations (`ChangeLog.md`, 2026-08-16,
    `SR.m`); this port follows v4.0. `Ball`, `IRall`, `VDall` and every
    band are unaffected — they depend only on the drawn `F` and `sigma`.

!!! warning "Throws when the draw budget is exhausted"
    If `VARopt.sr_draw` total draws pass without collecting `ndraws`
    accepted rotations, this raises. That means the restrictions are too
    tight, not that anything is broken — loosen `SIGN`, reduce `sr_hor`, or
    raise `sr_draw`.

!!! note "`SRout.IR`/`VD`/`HD` correspond to `SRout.B` (upstream v4.0 fix)"
    VAR Toolbox 3.1's `SR.m` stored each draw's model under a label built
    from `jj` **before** `jj` was incremented, so accepted draw `k` was
    filed under `"draw\$(k-1)"`, while the final lookup used `sel` directly:
    `SRout.IR`, `SRout.VD` and `SRout.HD` were computed from the accepted
    draw *after* the one `SRout.B` came from (verified in MATLAB on
    Stock-Watson data: `B` was draw 24, `IR` matched draw 25 exactly), and
    `sel == ndraws` crashed on a missing field. This port found and
    reported the defect; v4.0 `SR.m` now reads `"draw\$(sel-1)"`, and the
    default here matches v4.0. Pass `fix_draw_index = false` to reproduce
    the v3.1 mapping (used by the golden-master tests against v3.1
    fixtures).

The `med`/`inf`/`sup` bands are unaffected either way: they are computed
from `IRall`/`VDall` across all draws.
"""
function SR(VAR::VARModel, SIGN::AbstractMatrix, VARopt::VAROptions;
            rng::AbstractRNG = Random.default_rng(),
            src::Union{Nothing,DrawSource} = nothing,
            post_draws = nothing,
            fix_draw_index::Bool = true)

    src === nothing && (src = DrawSource(rng = rng))

    nvar, nvar_ex = VAR.nvar, VAR.nvar_ex
    nsteps, ndraws, pctg = VARopt.nsteps, VARopt.ndraws, VARopt.pctg
    nobs, nlag = VAR.nobs, VAR.nlag

    IRall = fill(NaN, nsteps, nvar, nvar, ndraws)
    VDall = fill(NaN, nsteps, nvar, nvar, ndraws)
    Ball = fill(NaN, nvar, nvar, ndraws)
    HDall = Dict{Symbol,Array{Float64}}(
        :shock => zeros(nobs + nlag, nvar, nvar, ndraws),
        :init => zeros(nobs + nlag, nvar, ndraws),
        :constant => zeros(nobs + nlag, nvar, ndraws),
        :trend => zeros(nobs + nlag, nvar, ndraws),
        :trend2 => zeros(nobs + nlag, nvar, ndraws),
        :endo => zeros(nobs + nlag, nvar, ndraws),
        :exo => zeros(nobs + nlag, nvar, max(nvar_ex, 1), ndraws),
    )

    # Keyed exactly as SR.m keys them: by the accepted count BEFORE the
    # increment. See the off-by-one note in the docstring.
    VAR_draw = Dict{Int,VARModel}()

    opt = @set VARopt.ident = "sign"

    jj = 0      # accepted draws
    tt = 0      # total draws
    npost = 0
    while jj < ndraws
        if tt > VARopt.sr_draw
            error("Total number of draws for finding sign restrictions exceeded " *
                  "($(VARopt.sr_draw)). Change the restrictions " *
                  "(SIGN / sr_hor) or increase VARopt.sr_draw.")
        end

        label = jj
        VARd = VAR
        if VARopt.sr_mod == 1
            npost += 1
            sd = post_draws === nothing ? nothing : post_draws[npost][1]
            fd = post_draws === nothing ? nothing : post_draws[npost][2]
            sigma_draw, Ft_draw, F_draw, Fcomp_draw =
                VARdrawpost(VAR; rng = rng, sigma_draw = sd, Ft_draw = fd)
            VARd = @set VARd.Ft = Ft_draw
            VARd = @set VARd.F = F_draw
            VARd = @set VARd.Fcomp = Fcomp_draw
            VARd = @set VARd.sigma = sigma_draw
            # v4.0 (SR.m): the drawn model carries its own innovations.
            # v3.1 left `resid` at the OLS point estimate.
            VARd = @set VARd.resid = VAR.Y - VAR.X * Ft_draw
        end

        B = SignRestrictions(SIGN, VARd, opt; src = src)

        if B !== nothing
            jj += 1
            tt += 1
            Ball[:, :, jj] = B
            VARd = @set VARd.B = B

            aux_irf, VARd = VARir(VARd, opt)
            IRall[:, :, :, jj] = aux_irf
            aux_fevd, _ = VARvd(VARd, opt)
            VDall[:, :, :, jj] = aux_fevd
            aux_hd, _ = VARhd(VARd, opt)
            HDall[:shock][:, :, :, jj] = aux_hd.shock
            HDall[:init][:, :, jj] = aux_hd.init
            HDall[:constant][:, :, jj] = aux_hd.constant
            HDall[:trend][:, :, jj] = aux_hd.trend
            HDall[:trend2][:, :, jj] = aux_hd.trend2
            HDall[:endo][:, :, jj] = aux_hd.endo
            nvar_ex > 0 && (HDall[:exo][:, :, :, jj] = aux_hd.exo)

            VAR_draw[label] = VARd
        else
            tt += 1
        end
    end

    Bmed = dropdims(median(Ball, dims = 3), dims = 3)
    aux = dropdims(sum((Ball .- Bmed) .^ 2, dims = (1, 2)), dims = (1, 2))
    sel = argmin(aux)

    pctg_inf = (100 - pctg) / 2
    pctg_sup = 100 - (100 - pctg) / 2

    IRmed = dropdims(median(IRall, dims = 4), dims = 4)
    IRb = _prctile(IRall, [pctg_inf, pctg_sup], 4)
    VDmed = dropdims(median(VDall, dims = 4), dims = 4)
    VDb = _prctile(VDall, [pctg_inf, pctg_sup], 4)

    # v4.0 SR.m reads `VAR_draw.('draw' + (sel-1))`, the draw B came from.
    # v3.1 read `'draw' + sel`, the draw AFTER it; `fix_draw_index=false`
    # reproduces that mapping for comparison against v3.1 fixtures.
    key = fix_draw_index ? sel - 1 : sel
    haskey(VAR_draw, key) || error(
        "SR: no stored draw for index $key. This reproduces the VAR Toolbox " *
        "3.1 defect in SR.m -- it labelled stored draws by the accepted count " *
        "before incrementing it, so the label for the last accepted draw " *
        "was never created and `sel == ndraws` read a missing field. " *
        "Use the default fix_draw_index=true.")
    VARsel = VAR_draw[key]

    IR, _ = VARir(VARsel, opt)
    VD, _ = VARvd(VARsel, opt)
    HD, _ = VARhd(VARsel, opt)

    return SROutput(
        IRall = IRall, VDall = VDall, Ball = Ball, HDall = HDall,
        Bmed = Bmed, B = Ball[:, :, sel],
        IRmed = IRmed,
        IRinf = IRb[:, :, :, 1], IRsup = IRb[:, :, :, 2],
        VDmed = VDmed,
        VDinf = VDb[:, :, :, 1], VDsup = VDb[:, :, :, 2],
        IR = IR, VD = VD, HD = HD,
    )
end
