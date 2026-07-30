# =======================================================================
# bands.jl
#
# Port of v3dot1/VAR/VARirband.m and VARvdband.m.
#
# Both re-estimate the VAR on artificial data over VARopt.ndraws bootstrap
# replications. Randomness goes through a DrawSource (M2) so MATLAB's own
# draws can be injected and the bands compared exactly.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

"""
    next_uniform!(src, dims...)

Uniform `[0,1)` draws from a [`DrawSource`](@ref).

Kept separate from `next_normal!` so an injected pool cannot be silently
consumed by the wrong kind of draw: the bootstrap uses `rand`, the
sign-restriction rotations use `randn`, and a pool recorded from one is
meaningless to the other.
"""
function next_uniform!(s::DrawSource, dims::Integer...)
    if s.uniforms === nothing
        return rand(s.rng, dims...)
    end
    s.j += 1
    s.j > length(s.uniforms) &&
        error("DrawSource: injected uniform pool exhausted after $(s.j - 1) draws")
    a = s.uniforms[s.j]
    size(a) == dims ||
        error("DrawSource: injected uniform $(s.j) has size $(size(a)), expected $(dims)")
    return copy(a)
end

# Generate one bootstrap sample of artificial data. Shared verbatim between
# VARirband.m and VARvdband.m (their artificial-data blocks are identical).
#
# Returns (y_artificial, Z) where Z is the resampled instrument, or nothing
# when the wild+iv branch was not taken.
function _artificial_data(VAR::VARModel, VARopt::VAROptions, src::DrawSource;
                          allow_iv_wild::Bool)
    method = VARopt.method
    nvar, nlag = VAR.nvar, VAR.nlag
    nobs = VAR.nobs
    resid = VAR.resid
    ENDO = VAR.ENDO
    constant = VAR.constant
    Ft = VAR.Ft

    Z = nothing
    if method == "bs"
        idx = ceil.(Int, size(resid, 1) .* vec(next_uniform!(src, nobs, 1)))
        u = resid[idx, :]
    elseif method == "wild"
        if allow_iv_wild && VARopt.ident == "iv"
            nIV = size(VAR.IV, 2)
            rr = 1 .- 2 .* (next_uniform!(src, nobs, nIV) .> 0.5)
            u = resid .* (rr * ones(nIV, nvar))
            Z = vcat(VAR.IV[1:nlag, :], VAR.IV[(nlag+1):end, :] .* rr)
        else
            rr = 1 .- 2 .* (next_uniform!(src, nobs, 1) .> 0.5)
            u = resid .* (rr * ones(1, nvar))
        end
    else
        throw(ArgumentError("The method $method is not available"))
    end

    y = zeros(nobs + nlag, nvar)
    LAG = Float64[]
    for jj in 1:nlag
        y[jj, :] = ENDO[jj, :]
        LAG = vcat(y[jj, :], LAG)
    end
    T = collect(1.0:nobs)

    function lagplus(LAG, t_index)
        if constant == 0
            v = copy(LAG)
        elseif constant == 1
            v = vcat(1.0, LAG)
        elseif constant == 2
            v = vcat(1.0, T[t_index], LAG)
        else
            v = vcat(1.0, T[t_index], T[t_index]^2, LAG)
        end
        VAR.nvar_ex != 0 && (v = vcat(v, VAR.X_EX[t_index, :]))
        return v
    end

    LAGplus = lagplus(LAG, 1)
    for jj in (nlag+1):(nobs+nlag)
        for mm in 1:nvar
            y[jj, mm] = dot(LAGplus, Ft[:, mm]) + u[jj-nlag, mm]
        end
        if jj < nobs + nlag
            LAG = vcat(y[jj, :], LAG[1:((nlag-1)*nvar)])
            LAGplus = lagplus(LAG, jj - nlag + 1)
        end
    end
    return y, Z
end

"""
    VARirband(VAR, VARopt; rng=Random.default_rng(), src=nothing)
        -> (INF, SUP, MED, BAR)

Bootstrap error bands for impulse responses. Each of the four outputs is
`(nsteps × nvar × nshocks)`, matching [`VARir`](@ref)'s layout.

`VARopt.method` selects `"bs"` (resample residuals with replacement) or
`"wild"` (Rademacher wild bootstrap). `VARopt.pctg` sets the coverage;
`VARopt.ndraws` the number of replications.

Draws whose companion matrix has `maxEig >= 0.9999` are **discarded and
redrawn**, so `ndraws` counts accepted replications, not attempts. There is
no cap on attempts — an explosive point estimate will loop indefinitely.
That is upstream behaviour.

`MED` is the pointwise median across draws and `BAR` the pointwise mean;
neither equals the point estimate from `VARir`.

!!! note "Percentiles use MATLAB's definition"
    Bands come from [`_prctile`](@ref), the Hazen/type-5 quantile MATLAB's
    `prctile` implements — not Julia's default. At the toolbox default
    `pctg = 95` the two differ by percent-level amounts. See M1.

Pass `src` (a [`DrawSource`](@ref) carrying recorded uniforms) to reproduce
a MATLAB run exactly.
"""
function VARirband(VAR::VARModel, VARopt::VAROptions;
                   rng::AbstractRNG = Random.default_rng(),
                   src::Union{Nothing,DrawSource} = nothing)
    src === nothing && (src = DrawSource(rng = rng))
    nsteps, ndraws, pctg = VARopt.nsteps, VARopt.ndraws, VARopt.pctg
    nvar = VAR.nvar

    IR = zeros(nsteps, nvar, nvar, ndraws)
    tt = 1
    while tt <= ndraws
        y, Z = _artificial_data(VAR, VARopt, src; allow_iv_wild = true)
        VAR_draw, _ = VAR.nvar_ex != 0 ?
                      VARmodel(y, VAR.nlag, VAR.constant, VAR.EXOG, VAR.nlag_ex) :
                      VARmodel(y, VAR.nlag, VAR.constant)
        Z !== nothing && (VAR_draw = @set VAR_draw.IV = Z)

        IR_draw, VAR_draw = VARir(VAR_draw, VARopt)
        if VAR_draw.maxEig < 0.9999
            IR[:, :, :, tt] = IR_draw
            tt += 1
        end
    end

    pctg_inf = (100 - pctg) / 2
    pctg_sup = 100 - (100 - pctg) / 2
    INF = dropdims(_prctile(IR, pctg_inf, 4), dims = 4)
    SUP = dropdims(_prctile(IR, pctg_sup, 4), dims = 4)
    MED = dropdims(_prctile(IR, 50, 4), dims = 4)
    BAR = dropdims(mean(IR, dims = 4), dims = 4)
    return INF, SUP, MED, BAR
end

"""
    VARvdband(VAR, VARopt; rng=Random.default_rng(), src=nothing)
        -> (INF, SUP, MED, BAR)

Bootstrap error bands for the forecast error variance decomposition. Each
output is `(nsteps × nshocks × nvar)`, matching [`VARvd`](@ref)'s layout —
which is **not** [`VARirband`](@ref)'s.

!!! warning "No IV wild-bootstrap branch, unlike `VARirband`"
    `VARirband.m` has a dedicated `wild` path for `ident = "iv"` that also
    resamples the instrument
    (`Z = [IV(1:nlag,:); IV(nlag+1:end,:).*rr]`). `VARvdband.m` has no such
    branch: with `method = "wild"` it always takes the plain
    `rand(nobs,1)` path regardless of `ident`.

    That asymmetry is replicated rather than fixed. In practice it rarely
    bites, because `VARvd` refuses `ident = "iv"` outright — but it means
    the two band functions consume different numbers of random draws for
    the same options, which matters when injecting recorded draws.

    (`VARvdband.m` also hardcodes its progress interval as `10*ww` instead
    of reading `VARopt.mult`. Cosmetic; not reproduced, since this port
    prints no progress.)
"""
function VARvdband(VAR::VARModel, VARopt::VAROptions;
                   rng::AbstractRNG = Random.default_rng(),
                   src::Union{Nothing,DrawSource} = nothing)
    src === nothing && (src = DrawSource(rng = rng))
    nsteps, ndraws, pctg = VARopt.nsteps, VARopt.ndraws, VARopt.pctg
    nvar = VAR.nvar

    VD = zeros(nsteps, nvar, nvar, ndraws)
    tt = 1
    while tt <= ndraws
        y, _ = _artificial_data(VAR, VARopt, src; allow_iv_wild = false)
        VAR_draw, _ = VAR.nvar_ex != 0 ?
                      VARmodel(y, VAR.nlag, VAR.constant, VAR.EXOG, VAR.nlag_ex) :
                      VARmodel(y, VAR.nlag, VAR.constant)

        vd_draw, VAR_draw_opt = VARvd(VAR_draw, VARopt)
        if VAR_draw_opt.maxEig < 0.9999
            VD[:, :, :, tt] = vd_draw
            tt += 1
        end
    end

    pctg_inf = (100 - pctg) / 2
    pctg_sup = 100 - (100 - pctg) / 2
    INF = dropdims(_prctile(VD, pctg_inf, 4), dims = 4)
    SUP = dropdims(_prctile(VD, pctg_sup, 4), dims = 4)
    MED = dropdims(_prctile(VD, 50, 4), dims = 4)
    BAR = dropdims(mean(VD, dims = 4), dims = 4)
    return INF, SUP, MED, BAR
end
