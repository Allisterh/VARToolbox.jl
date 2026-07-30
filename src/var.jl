# =======================================================================
# var.jl
#
# Port of v3dot1/VAR/VARoption.m and v3dot1/VAR/VARmodel.m.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

# -----------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------

"""
    VAROptions

Flat options struct controlling everything downstream of [`VARmodel`](@ref):
identification, error bands, and (once the plotting extension lands) figures.
Constructed with defaults by [`VARoption`](@ref).

!!! important "Immutable by design (D3)"
    MATLAB passes structs **by value**, and the toolbox depends on it:
    `SR.m` sets `VARopt.ident = 'sign'` and `SignRestrictions.m` sets
    `VARopt.nsteps = VARopt.sr_hor` on their local copies, expecting the
    caller's options to be untouched. A mutable Julia struct would leak
    those mutations back and silently corrupt every band and plot computed
    afterwards.

    Update by reconstruction with `@set` (re-exported by this package):

    ```julia
    VARopt = VARoption()
    VARopt = @set VARopt.nsteps = 24
    VARopt = @set VARopt.ident  = "short"
    ```

Fields mirror `VARoption.m` one-for-one.

| field | default | meaning |
|---|---|---|
| `vnames`, `vnames_ex`, `snames` | `String[]` | endogenous / exogenous / shock names |
| `nsteps` | `40` | horizon for IRFs and FEVDs |
| `impact` | `0` | shock size: `0` one s.d., `1` unit |
| `shut` | `0` | force the IRF of one variable to zero |
| `ident` | `"short"` | `"short"`, `"long"`, `"sign"`, `"iv"`, `"exog"` |
| `recurs` | `"wold"` | `"wold"` MA representation, `"comp"` companion form |
| `ndraws` | `1000` | bootstrap / sign-restriction draws |
| `mult` | `10` | progress-printing interval |
| `pctg` | `95` | confidence level |
| `method` | `"bs"` | `"bs"` standard bootstrap, `"wild"` wild bootstrap |
| `sr_hor` | `1` | horizons over which sign restrictions bind |
| `sr_rot` | `500` | max rotations per draw |
| `sr_draw` | `100000` | max total draws |
| `sr_mod` | `1` | include model uncertainty (redraw F, sigma) |
| `pick` | `0` | plot a single variable (`0` plots all) |
| `quality` | `2` | figure export quality |
| `suptitle` | `0` | title above figures |
| `datesnum`, `datestxt`, `datestype` | `Float64[]`, `String[]`, `1` | date labelling |
| `firstdate` | `nothing` | e.g. `1999.75` for 1999Q4 |
| `frequency` | `"q"` | `"m"`, `"q"`, `"y"` |
| `figname` | `""` | exported figure name |
| `FigSize` | `[26, 24]` | figure window size |
"""
Base.@kwdef struct VAROptions
    vnames::Vector{String} = String[]
    vnames_ex::Vector{String} = String[]
    snames::Vector{String} = String[]
    nsteps::Int = 40
    impact::Int = 0
    shut::Int = 0
    ident::String = "short"
    recurs::String = "wold"
    ndraws::Int = 1000
    mult::Int = 10
    pctg::Float64 = 95.0
    method::String = "bs"
    sr_hor::Int = 1
    sr_rot::Int = 500
    sr_draw::Int = 100000
    sr_mod::Int = 1
    pick::Int = 0
    quality::Int = 2
    suptitle::Int = 0
    datesnum::Vector{Float64} = Float64[]
    datestxt::Vector{String} = String[]
    datestype::Int = 1
    firstdate::Union{Nothing,Float64} = nothing
    frequency::String = "q"
    figname::String = ""
    FigSize::Vector{Int} = [26, 24]
end

"""
    VARoption() -> VAROptions

Default options, matching `VARoption.m`. Called automatically by
[`VARmodel`](@ref), which returns the result as its second output.
"""
VARoption() = VAROptions()

# Julia's default `==` for an immutable struct falls back to `===` on each
# field, so two independently-allocated empty `String[]` fields compare
# unequal. Compare by value instead -- users reasonably expect
# `VARoption() == VARoption()`.
function Base.:(==)(a::VAROptions, b::VAROptions)
    return all(f -> getfield(a, f) == getfield(b, f), fieldnames(VAROptions))
end
Base.hash(o::VAROptions, h::UInt) =
    foldr((f, acc) -> hash(getfield(o, f), acc), fieldnames(VAROptions); init = h)

# -----------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------

"""
    EquationResult

Per-equation OLS output from [`VARmodel`](@ref), reached as `VAR.eq[j]`.

Replaces the MATLAB pattern at `VARmodel.m:107-121`, which builds fields
named `VAR.eq1`, `VAR.eq2`, … through 11 `eval` calls, and is read back by
`VARprint.m` through 5 more. A typed vector removes all 16 while keeping the
access pattern recognisable.

!!! note "`tprob` uses different degrees of freedom from `OLSmodel`"
    `VARmodel` recomputes `tprob` as `tdis_prb(tstat, nobs - ncoeff)`, where
    `ncoeff = nvar*nlag` — i.e. it counts only the lag coefficients,
    excluding the constant and any exogenous regressors. `OLSmodel` reports
    `tdis_prb(tstat, nobs)`. Both are upstream behaviour; they disagree, and
    the `VARmodel` value is the one carried here.
"""
struct EquationResult
    beta::Vector{Float64}
    tstat::Vector{Float64}
    bstd::Vector{Float64}
    tprob::Vector{Float64}
    resid::Vector{Float64}
    yhat::Vector{Float64}
    y::Vector{Float64}
    rsqr::Float64
    rbar::Float64
    sige::Float64
    dw::Float64
end

"""
    VARModel

Estimation output from [`VARmodel`](@ref).

Immutable, for the same reason as [`VAROptions`](@ref): in MATLAB,
`VARir`/`VARvd`/`VARhd` receive a *copy*, fill in `B`/`PSI`/`Fp`, and return
it. The caller's original stays empty unless it captures the second return
value. Reproducing that with an immutable struct makes the contract explicit
rather than accidental:

```julia
IR, VAR = VARir(VAR, VARopt)   # capture the updated model
```

| field | shape | meaning |
|---|---|---|
| `ENDO`, `EXOG` | `T x nvar`, `T x nvar_ex` | input data |
| `nlag`, `nlag_ex`, `constant` | | lag orders and deterministics (`const` renamed) |
| `nobs` | | **effective** observations, `T - max(nlag, nlag_ex)` |
| `nvar`, `nvar_ex` | | variable counts |
| `ncoeff` | | `nvar*nlag` — lag coefficients only |
| `ntotcoeff` | | `ncoeff + nvar_ex*(nlag_ex+1) + constant` |
| `eq` | `nvar`-vector | per-equation results, see [`EquationResult`](@ref) |
| `Ft`, `F` | `k x nvar`, `nvar x k` | coefficients; `F == Ft'` |
| `sigma` | `nvar x nvar` | residual covariance, `nobs - ntotcoeff` denominator |
| `resid` | `nobs x nvar` | residuals |
| `X`, `Y`, `X_EX` | | design matrices |
| `Fcomp` | | companion matrix |
| `maxEig` | | largest absolute eigenvalue of `Fcomp` — stability check |
| `B`, `Biv` | | structural impact matrix; empty until identified |
| `PSI`, `Fp` | | Wold multipliers, recursive `F` by lag; empty until `VARir` |
| `IV` | | external instruments; empty until set |
"""
Base.@kwdef struct VARModel
    ENDO::Matrix{Float64}
    EXOG::Matrix{Float64} = zeros(0, 0)
    nlag::Int
    nlag_ex::Int = 0
    constant::Int = 1
    nobs::Int
    nvar::Int
    nvar_ex::Int = 0
    ncoeff::Int
    ntotcoeff::Int
    eq::Vector{EquationResult}
    Ft::Matrix{Float64}
    F::Matrix{Float64}
    sigma::Matrix{Float64}
    resid::Matrix{Float64}
    X::Matrix{Float64}
    Y::Matrix{Float64}
    X_EX::Matrix{Float64} = zeros(0, 0)
    Fcomp::Matrix{Float64}
    maxEig::Float64
    B::Matrix{Float64} = zeros(0, 0)
    Biv::Matrix{Float64} = zeros(0, 0)
    PSI::Array{Float64,3} = zeros(0, 0, 0)
    Fp::Array{Float64,3} = zeros(0, 0, 0)
    IV::Matrix{Float64} = zeros(0, 0)
    # Populated only by `ident = "iv"`, mirroring the extra fields VARir.m
    # attaches to VAR in that branch.
    FirstStage::Union{Nothing,OLSOutput} = nothing
    sigma_b::Matrix{Float64} = zeros(0, 0)
end

# -----------------------------------------------------------------------
# Estimation
# -----------------------------------------------------------------------

"""
    VARmodel(ENDO, nlag, const_=1, EXOG=nothing, nlag_ex=0) -> (VAR, VARopt)

Estimate a reduced-form VAR by OLS, equation by equation.

Returns the fitted [`VARModel`](@ref) **and** a default [`VAROptions`](@ref).
Everything downstream is configured by updating that options struct with
`@set` between calls — there are no keyword arguments anywhere in the
pipeline, matching the MATLAB API.

`const_` selects deterministics: `0` none, `1` constant, `2` constant plus
linear trend, `3` constant plus linear and squared trend. Named `const_`
because `const` is reserved in Julia.

!!! note "`const_` defaults to 1"
    `VARmodel.m`'s own docstring says `[dflt = 0]`, but the code defaults to
    `1` (constant included). The code is authoritative; the documentation
    bug is fixed here rather than reproduced — this is one of the two
    agreed carve-outs from the replicate-everything policy.

!!! warning "No identification is performed here"
    `B`, `Biv`, `PSI`, `Fp` and `IV` come back empty. The identification
    scheme is late-bound: set `VARopt.ident` and call [`VARir`](@ref),
    `VARvd` or `VARhd`, each of which returns an updated model as its second
    output. For sign restrictions go through `SR` instead.

# Example
```julia
VAR, VARopt = VARmodel(ENDO, 4)
VARopt = @set VARopt.nsteps = 24
IR, VAR = VARir(VAR, VARopt)
```
"""
function VARmodel(ENDO::AbstractMatrix, nlag::Integer, const_::Integer = 1,
                  EXOG::Union{Nothing,AbstractMatrix} = nothing, nlag_ex::Integer = 0)
    nlag < 1 && throw(ArgumentError("VARmodel: nlag must be >= 1 (got $nlag)"))
    nobs, nvar = size(ENDO)

    if EXOG !== nothing
        size(EXOG, 1) == nobs ||
            throw(DimensionMismatch("var: nobs in EXOG-matrix not the same as y-matrix"))
        nvar_ex = size(EXOG, 2)
    else
        nvar_ex = 0
        nlag_ex = 0
    end

    nobse = nobs - max(nlag, nlag_ex)
    ncoeff = nvar * nlag
    ncoeff_ex = nvar_ex * (nlag_ex + 1)
    ntotcoeff = ncoeff + ncoeff_ex + const_

    Y, X = VARmakexy(Matrix{Float64}(ENDO), nlag, const_)

    X_EX = zeros(0, 0)
    if nvar_ex > 0
        X_EX = VARmakelags(Matrix{Float64}(EXOG), nlag_ex)
        if nlag == nlag_ex
            X = hcat(X, X_EX)
        elseif nlag > nlag_ex
            d = nlag - nlag_ex
            X_EX = X_EX[(d+1):end, :]
            X = hcat(X, X_EX)
        else
            d = nlag_ex - nlag
            Y = Y[(d+1):end, :]
            X = hcat(X[(d+1):end, :], X_EX)
        end
    end

    # ---- per-equation OLS (replaces the 11 eval writes at VARmodel.m:107) --
    eqs = Vector{EquationResult}(undef, nvar)
    for j in 1:nvar
        yv = Y[:, j]
        O = OLSmodel(yv, X, 0)
        # Recomputed with nobse - ncoeff dof, NOT OLSmodel's nobs. Upstream.
        tprob = tdis_prb(O.tstat, nobse - ncoeff)
        eqs[j] = EquationResult(O.beta, O.tstat, O.bstd, tprob, O.resid,
                                O.yhat, yv, O.rsqr, O.rbar, O.sige, O.dw)
    end

    # ---- system coefficients -------------------------------------------
    # MATLAB's `(X'*X)\(X'*Y)` dispatches on the symmetry of X'X and solves
    # by Cholesky; Julia's generic `\` on a Matrix uses LU with partial
    # pivoting. Both are correct, but they round differently. Matching the
    # algorithm keeps agreement with the fixtures near machine precision
    # rather than merely within tolerance.
    XtX = transpose(X) * X
    XtY = transpose(X) * Y
    Ft = cholesky(Symmetric(XtX)) \ XtY
    F = Matrix(transpose(Ft))
    resid = Y - X * Ft
    sigma = (1 / (nobse - ntotcoeff)) .* (transpose(resid) * resid)

    # ---- companion form ------------------------------------------------
    k = nvar * (nlag - 1)
    top = F[:, (1+const_):(nvar*nlag+const_)]
    Fcomp = k == 0 ? top : vcat(top, hcat(Matrix{Float64}(I, k, k), zeros(k, nvar)))
    maxEig = maximum(abs.(eigvals(Fcomp)))

    VAR = VARModel(
        ENDO = Matrix{Float64}(ENDO),
        EXOG = EXOG === nothing ? zeros(0, 0) : Matrix{Float64}(EXOG),
        nlag = nlag, nlag_ex = nlag_ex, constant = const_,
        nobs = nobse, nvar = nvar, nvar_ex = nvar_ex,
        ncoeff = ncoeff, ntotcoeff = ntotcoeff,
        eq = eqs, Ft = Ft, F = F, sigma = sigma, resid = resid,
        X = X, Y = Y, X_EX = X_EX, Fcomp = Fcomp, maxEig = maxEig,
    )
    return VAR, VARoption()
end
