# =======================================================================
# ardl.jl
#
# Port of v3dot1/VAR/ARDLmodel.m -- a third, independent single-equation
# model class alongside VARmodel and LPmodel. It shares VARmodel's argument
# shape but is not part of the VAR/VARopt pipeline.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

"""
    ARDLModel

Result of [`ARDLmodel`](@ref). Field names follow `ARDLmodel.m`, with
`const` renamed to `constant` (reserved word).

!!! note "`nvar` here counts coefficients, not variables"
    Upstream sets `ARDL.nvar = ncoeff + ncoeff_ex + const`, i.e. the total
    number of estimated coefficients. The name is misleading but is kept for
    API parity. `nvar_ex` *does* count exogenous variables.
"""
Base.@kwdef struct ARDLModel
    meth::String = "ols"
    ENDO::Matrix{Float64}
    EXOG::Matrix{Float64} = zeros(0, 0)
    nlag::Int
    nlag_ex::Int = 0
    constant::Int = 1
    nobs::Int
    nvar::Int
    nvar_ex::Int = 0
    ncoeff::Int
    y::Matrix{Float64}
    x::Matrix{Float64}
    beta::Matrix{Float64}
    yhat::Matrix{Float64}
    resid::Matrix{Float64}
    sige::Float64
    sigbeta::Matrix{Float64}
    bstd::Vector{Float64}
    bint::Matrix{Float64}
    tstat::Vector{Float64}
    tprob::Vector{Float64}
    rsqr::Float64
    rbar::Float64
    dw::Float64
    F::Float64 = NaN
end

"""
    ARDLmodel(ENDO, nlag, const_=1, EXOG=nothing, nlag_ex=0) -> ARDLModel

Single-equation autoregressive distributed-lag model. `ENDO` must have one
column.

Same argument shape as [`VARmodel`](@ref), but this is a standalone class:
it returns no options struct and does not feed `VARir`/`VARvd`/`VARhd`.

!!! note "Degrees of freedom (upstream v4.0 fix)"
    `bint` and `tprob` both use the residual degrees of freedom
    `nobse - nvar` (effective sample less coefficients), the denominator of
    `sige`. VAR Toolbox 3.1 used `nobse` for the critical value and the
    *raw* `nobs` for the p-values — two conventions that disagreed with each
    other and with the residual variance; v4.0 `ARDLmodel.m` uses
    `nobse - nvar` for both. This port follows v4.0.
"""
function ARDLmodel(ENDO::AbstractMatrix, nlag::Integer, const_::Integer = 1,
                   EXOG::Union{Nothing,AbstractMatrix} = nothing, nlag_ex::Integer = 0)
    size(ENDO, 2) == 1 ||
        throw(ArgumentError("ARDLmodel: ENDO must have exactly one column (got $(size(ENDO,2)))"))
    nobs = size(ENDO, 1)

    if EXOG !== nothing
        size(EXOG, 1) == nobs ||
            throw(DimensionMismatch("var: nobs in EXOG-matrix not the same as y-matrix"))
        nvar_ex = size(EXOG, 2)
    else
        nvar_ex = 0
        nlag_ex = 0
    end

    nobse = nobs - max(nlag, nlag_ex)
    ncoeff = nlag
    ncoeff_ex = nvar_ex + nvar_ex * nlag_ex
    nvar = ncoeff + ncoeff_ex + const_

    Y, X = VARmakexy(Matrix{Float64}(ENDO), nlag, const_)

    if EXOG !== nothing
        X_EX = VARmakelags(Matrix{Float64}(EXOG), nlag_ex)
        if nlag == nlag_ex
            X = hcat(X, X_EX)
        elseif nlag > nlag_ex
            d = nlag - nlag_ex
            X = hcat(X, X_EX[(d+1):end, :])
        else
            d = nlag_ex - nlag
            Y = Y[(d+1):end, :]
            X = hcat(X[(d+1):end, :], X_EX)
        end
    end

    R = qr(X).R
    xpxi = inv(transpose(R) * R)

    beta = xpxi * (transpose(X) * Y)
    yhat = X * beta
    resid = Y - yhat

    sigu = (transpose(resid)*resid)[1, 1]
    sige = sigu / (nobse - nvar)
    sigbeta = sige .* xpxi

    tmp = sige .* diag(xpxi)
    bstd = sqrt.(tmp)
    # Residual degrees of freedom (v4.0); v3.1 used nobse and nobs here.
    tcrit = -tdis_inv(0.025, nobse - nvar)
    bint = hcat(vec(beta) .- tcrit .* bstd, vec(beta) .+ tcrit .* bstd)
    tstat = vec(beta) ./ sqrt.(tmp)
    tprob = tdis_prb(tstat, nobse - nvar)

    ym = Y .- mean(Y)
    rsqr2raw = (transpose(ym)*ym)[1, 1]
    rsqr = 1.0 - sigu / rsqr2raw
    rsqr1 = sigu / (nobse - nvar)
    rsqr2 = rsqr2raw / (nobse - 1.0)
    rbar = rsqr2 != 0 ? 1 - (rsqr1 / rsqr2) : rsqr

    ediff = resid[2:nobse, 1] - resid[1:(nobse-1), 1]
    dw = dot(ediff, ediff) / sigu

    Fstat = NaN
    if const_ > 0
        fx = X[:, 1:1]
        fbeta = (transpose(fx) * fx) \ (transpose(fx) * Y)
        fresid = Y - fx * fbeta
        frsqr = 1.0 - (transpose(fresid)*fresid)[1, 1] / rsqr2raw
        Fstat = ((frsqr - rsqr) / (1 - nvar)) / ((1 - rsqr) / (nobse - nvar))
    end

    return ARDLModel(
        ENDO = Matrix{Float64}(ENDO),
        EXOG = EXOG === nothing ? zeros(0, 0) : Matrix{Float64}(EXOG),
        nlag = nlag, nlag_ex = nlag_ex, constant = const_,
        nobs = nobse, nvar = nvar, nvar_ex = nvar_ex, ncoeff = ncoeff,
        y = Y, x = X, beta = beta, yhat = yhat, resid = resid,
        sige = sige, sigbeta = sigbeta, bstd = bstd, bint = bint,
        tstat = tstat, tprob = tprob, rsqr = rsqr, rbar = rbar, dw = dw, F = Fstat,
    )
end
