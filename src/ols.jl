# =======================================================================
# ols.jl
#
# Port of v3dot1/VAR/OLSmodel.m -- the single-equation OLS engine shared by
# VARmodel (called once per equation) and ARDLmodel.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

"""
    OLSOutput

Result of [`OLSmodel`](@ref). Field names follow `OLSmodel.m` with one
forced rename: MATLAB's `OLS.const` becomes `constant`, because `const` is
a reserved word in Julia.

| field | meaning |
|---|---|
| `meth` | `"ols"` |
| `y`, `x` | regressand, and regressors **including** any added deterministics |
| `nobs`, `nvar` | observations, and regressors after deterministics are added |
| `beta`, `yhat`, `resid` | coefficients, fit, residuals |
| `sige` | residual variance, `n - k` denominator |
| `sigbeta` | coefficient covariance |
| `bstd`, `tstat`, `tprob` | standard errors, t-statistics, two-sided p-values |
| `bint` | 95% confidence interval, two columns |
| `bstd_HW`, `bstd_NW` | White and Newey-West standard errors |
| `rsqr`, `rbar`, `dw` | R-squared, adjusted R-squared, Durbin-Watson |
| `constant` | the `const_` argument that was used |
| `F` | F-statistic; `NaN` when `const_ == 0` (upstream leaves it unset) |
"""
struct OLSOutput
    meth::String
    y::Vector{Float64}
    x::Matrix{Float64}
    nobs::Int
    nvar::Int
    beta::Vector{Float64}
    yhat::Vector{Float64}
    resid::Vector{Float64}
    sige::Float64
    sigbeta::Matrix{Float64}
    bstd::Vector{Float64}
    bint::Matrix{Float64}
    tstat::Vector{Float64}
    tprob::Vector{Float64}
    bstd_HW::Vector{Float64}
    bstd_NW::Vector{Float64}
    rsqr::Float64
    rbar::Float64
    dw::Float64
    constant::Int
    F::Float64
end

# Heteroskedasticity- and autocorrelation-consistent covariance, Bartlett
# weights. `nlag = 0` gives White; `nlag = floor(4*(n/100)^(2/9))` gives the
# Newey-West bandwidth upstream uses.
function _hac_std(resid::AbstractVector, x::AbstractMatrix, xpxi::AbstractMatrix,
                  nobs::Integer, nvar::Integer, nlag::Integer)
    hhat = repeat(reshape(resid, 1, :), nvar, 1) .* transpose(x)
    G = zeros(nvar, nvar)
    for a in 0:nlag
        w = (nlag + 1 - a) / (nlag + 1)
        za = @views hhat[:, (a+1):nobs] * transpose(hhat[:, 1:(nobs-a)])
        ga = a == 0 ? za : za + transpose(za)
        G .+= w .* ga
    end
    V = xpxi * G * xpxi
    return sqrt.(diag(V))
end

"""
    OLSmodel(y, x, const_=1) -> OLSOutput

Single-equation OLS.

`const_` selects deterministics: `0` none, `1` constant, `2` constant plus
linear trend.

!!! danger "`const_ == 3` is broken upstream and is rejected here"
    `OLSmodel.m`'s `const==3` branch prepends only two columns
    (`[ones trend.^2]`) while incrementing `nvar` by three, so the
    subsequent `xpxi = (r'*r)\\eye(nvar)` is built at the wrong dimension
    and MATLAB errors out. Note also that `VARmakexy.m`'s `const==3`
    prepends *three* columns (constant, trend, trend squared) -- the two
    functions disagree about what `3` means. Rather than replicate a
    guaranteed dimension-mismatch crash, this port throws an explicit
    error naming the upstream defect.

!!! note "Degrees of freedom in `tprob`"
    `tstat`'s p-values and the `bint` critical value both use `nobs`
    degrees of freedom, not `nobs - nvar`. That is upstream's choice and is
    replicated. [`VARmodel`](@ref) recomputes its own `tprob` with the
    correct `nobs - ncoeff`, so the per-equation `tprob` inside a fitted VAR
    differs from what `OLSmodel` alone reports.
"""
function OLSmodel(y::AbstractVector, x::AbstractMatrix, const_::Integer = 1)
    nobs = length(y)
    if !isempty(x)
        size(x, 1) == nobs || throw(DimensionMismatch("x and y must have same # obs"))
    end
    nvar = size(x, 2)

    X = Matrix{Float64}(x)
    if const_ == 1
        X = hcat(ones(nobs), X)
        nvar += 1
    elseif const_ == 2
        X = hcat(ones(nobs), collect(1.0:nobs), X)
        nvar += 2
    elseif const_ == 3
        throw(ArgumentError(
            "OLSmodel: const_ == 3 is broken in the MATLAB original " *
            "(OLSmodel.m adds 2 columns but increments nvar by 3, so the " *
            "normal-equations inverse is built at the wrong dimension). " *
            "Note VARmakexy's const==3 means constant + trend + trend^2. " *
            "Use const_ = 2, or build the regressors explicitly."))
    elseif const_ != 0
        throw(ArgumentError("OLSmodel: const_ must be 0, 1 or 2 (got $const_)"))
    end

    yv = Vector{Float64}(y)

    # Upstream switches to the normal equations above 10,000 observations;
    # both branches compute inv(X'X), so the QR path is used unconditionally
    # here for its better conditioning.
    R = qr(X).R
    xpxi = inv(transpose(R) * R)

    beta = xpxi * (transpose(X) * yv)
    yhat = X * beta
    resid = yv - yhat

    sigu = dot(resid, resid)
    sige = sigu / (nobs - nvar)
    sigbeta = sige .* xpxi

    tmp = sige .* diag(xpxi)
    sigb = sqrt.(tmp)
    tcrit = -tdis_inv(0.025, nobs)
    bint = hcat(beta .- tcrit .* sigb, beta .+ tcrit .* sigb)
    tstat = beta ./ sqrt.(tmp)
    tprob = tdis_prb(tstat, nobs)

    bstd_HW = _hac_std(resid, X, xpxi, nobs, nvar, 0)
    bstd_NW = _hac_std(resid, X, xpxi, nobs, nvar, floor(Int, 4 * ((nobs / 100)^(2 / 9))))

    ym = yv .- mean(yv)
    rsqr2raw = dot(ym, ym)
    rsqr = 1.0 - sigu / rsqr2raw
    rsqr1 = sigu / (nobs - nvar)
    rsqr2 = rsqr2raw / (nobs - 1.0)
    rbar = rsqr2 != 0 ? 1 - (rsqr1 / rsqr2) : rsqr

    ediff = resid[2:nobs] - resid[1:(nobs-1)]
    dw = dot(ediff, ediff) / sigu

    Fstat = NaN
    if const_ > 0
        fx = X[:, 1:1]
        fbeta = (transpose(fx) * fx) \ (transpose(fx) * yv)
        fresid = yv - fx * fbeta
        frsqr = 1.0 - dot(fresid, fresid) / rsqr2raw
        Fstat = ((frsqr - rsqr) / (1 - nvar)) / ((1 - rsqr) / (nobs - nvar))
    end

    return OLSOutput("ols", yv, X, nobs, nvar, beta, yhat, resid, sige, sigbeta,
                     sigb, bint, tstat, tprob, bstd_HW, bstd_NW,
                     rsqr, rbar, dw, const_, Fstat)
end
