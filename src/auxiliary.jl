# =======================================================================
# auxiliary.jl
#
# Port of v3dot1/Auxiliary/ -- the third-party and legacy numerical
# routines the toolbox depends on.
#
# Two things in here are load-bearing for numerical agreement with MATLAB
# and are documented at length below:
#
#   * `_prctile` implements MATLAB's `prctile`, which is the Hazen (type-5)
#     quantile definition, NOT Julia's default (type 7). See M1.
#   * `tdis_prb` / `tdis_inv` / `beta_inv` are hand-rolled upstream on base
#     MATLAB's `betainc`, and their clamping/iteration details perturb
#     results well above 1e-10. They are replicated, not substituted. See M3.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

# -----------------------------------------------------------------------
# prctile
# -----------------------------------------------------------------------

"""
    _prctile(x, p)
    _prctile(x, p, dim)

MATLAB-compatible `prctile`. `p` is in **percent** (0-100) and may be a
scalar or a vector.

!!! warning "This is not Julia's default quantile"
    MATLAB's `prctile` uses the Hazen / type-5 plotting position, i.e.
    `quantile(v, q; alpha=0.5, beta=0.5)`. Julia's `Statistics.quantile`
    defaults to type 7. On the toolbox's own default confidence level
    (`pctg = 95`) the two differ by percent-level amounts, which is far
    larger than any tolerance a bootstrap band test should use. Every
    `prctile` call site in the MATLAB source routes through this function.

NaNs are omitted, matching MATLAB. An all-NaN slice yields `NaN`.
"""
function _prctile end

function _prctile1(v::AbstractArray, p::Real)
    w = Float64[x for x in v if !isnan(x)]
    isempty(w) && return NaN
    return quantile(w, p / 100; alpha = 0.5, beta = 0.5)
end

_prctile1(v::AbstractArray, p::AbstractVector) = [_prctile1(v, pi) for pi in p]

_prctile(x::AbstractArray, p) = _prctile1(x, p)

function _prctile(x::AbstractArray, p, dim::Integer)
    dim > ndims(x) && return copy(x)
    return mapslices(v -> _prctile1(v, p), x; dims = dim)
end

# -----------------------------------------------------------------------
# Distribution functions
#
# Upstream these are hand-rolled on base MATLAB's `betainc` -- they do NOT
# require the Statistics and Machine Learning Toolbox. `betainc(x,a,b)` is
# the regularized incomplete beta, i.e. `cdf(Beta(a,b), x)`.
# -----------------------------------------------------------------------

_betainc(x, a, b) = cdf(Beta(a, b), x)

"""
    tdis_prb(x, n)

Two-sided p-value for a t-statistic `x` with `n` degrees of freedom, i.e.
`P(|T| > |x|)` for `T ~ t(n)`.

Analytically this equals `2*ccdf(TDist(n), abs(x))`, but the MATLAB source
clamps its beta argument to `1 - 1e-12` before calling `betainc`. Because
the `Beta(n/2, 1/2)` density is singular at 1, that clamp perturbs the
result by O(1e-6) near `x = 0` -- four orders of magnitude above the
tolerance used for golden-master tests. **The clamp is replicated
deliberately**; do not "simplify" this to the closed form.

The `x2 <= 0` clamp is unreachable for real `x` (since `x2 = n/(n+x^2) > 0`)
but is kept for fidelity.
"""
function tdis_prb(x, n::Real)
    n <= 0 && throw(ArgumentError("dof is negative or zero in tdis_prb"))
    return map(x) do xi
        x2 = n / (n + xi^2)
        x2 >= 1 && (x2 = 1 - 1e-12)     # replicated MATLAB clamp -- see docstring
        x2 <= 0 && (x2 = 1e-12)
        tmp = 1.0 - 0.5 * _betainc(x2, 0.5 * n, 0.5)
        2 * (1 - tmp)
    end
end

"""
    beta_pdf(x, a, b)

Beta density, zero outside `[0,1]`. Matches `Auxiliary/beta_pdf.m`.
"""
function beta_pdf(x, a::Real, b::Real)
    (a <= 0 || b <= 0) && throw(ArgumentError("Parameter a or b is nonpositive"))
    d = Beta(a, b)
    return map(xi -> (xi < 0 || xi > 1) ? zero(float(xi)) : pdf(d, xi), x)
end

"""
    beta_inv(p, a, b)

Inverse Beta CDF by the safeguarded Newton iteration in
`Auxiliary/beta_inv.m`.

Replicated rather than replaced by `quantile(Beta(a,b), p)`: the upstream
convergence test is `abs(dx) > 256*eps*max(x,1)`, so the iterate stops at a
point that differs from the exact quantile in the last few ulp. Replicating
keeps `tdis_inv` -- and hence `OLSmodel`'s confidence intervals -- bit-comparable
with MATLAB.
"""
function beta_inv(p, a::Real, b::Real)
    (a <= 0 || b <= 0) && throw(ArgumentError("beta_inv parameter a or b is nonpositive"))
    any(pi -> abs(2 * pi - 1) > 1, p) &&
        throw(ArgumentError("beta_inv: A probability should be 0<=p<=1"))

    x = fill(a / (a + b), size(p))
    pv = collect(float.(p))
    dx = fill(1.0, size(x))
    iter = 0
    while any(i -> abs(dx[i]) > 256 * eps() * max(x[i], 1), eachindex(x))
        iter += 1
        iter > 10_000 && break              # not upstream; guards a pathological input
        for i in eachindex(x)
            dx[i] = (_betainc(x[i], a, b) - pv[i]) / beta_pdf(x[i], a, b)
            x[i] = x[i] - dx[i]
            x[i] < 0 && (x[i] = x[i] + (dx[i] - x[i]) / 2)
        end
    end
    return p isa Number ? x[1] : x
end

"""
    tdis_inv(p, a)

Inverse t CDF with `a` degrees of freedom, via `beta_inv`. Matches
`Auxiliary/tdis_inv.m`.

!!! warning "Accuracy collapses near p = 0.5"
    The final step is `sqrt(x*a/(1-x))`, and `beta_inv` stops iterating once
    `|dx| <= 256*eps*max(x,1)`, i.e. at an absolute residual around 1e-14.
    Taking a square root of that residual amplifies it to ~1e-7. Concretely,
    `tdis_inv(0.5, 10)` returns `3.2148354647791181e-07` in MATLAB where the
    exact answer is `0`; this port reproduces that to 16 significant digits.

    The error is negligible in the tails, so `OLSmodel`'s use of
    `tdis_inv(0.025, nobs)` for confidence intervals is unaffected. If you
    need an accurate central quantile, use `quantile(TDist(a), p)`.
"""
function tdis_inv(p, a::Real)
    pv = float.(p)
    s = pv .< 0.5
    q = pv .+ (1 .- 2 .* pv) .* s          # reflect the lower tail
    q = 1 .- (2 .* (1 .- q))               # == 2q - 1
    x = beta_inv(q, 1 / 2, a / 2)
    x = x .* a ./ (1 .- x)
    return (1 .- 2 .* s) .* sqrt.(x)
end

"""
    fdis_prb(x, m, n)

Upper-tail probability of the F(m,n) distribution. Matches
`Auxiliary/fdis_prb.m`, including its `tst == 0` / `tst == 1` nudges of
1e-4 (which are far coarser than `tdis_prb`'s 1e-12 clamp -- upstream
inconsistency, replicated).
"""
function fdis_prb(x, m::Real, n::Real)
    (m <= 0 || n <= 0) && throw(ArgumentError("negative or zero dof in fdis_prb"))
    return map(x) do xi
        xi < 0 && throw(ArgumentError("negative or zero dof in fdis_prb"))
        tst = n / (n + m * xi)
        tst == 0 && (tst += 0.0001)
        tst == 1 && (tst -= 0.0001)
        _betainc(tst, 0.5 * n, 0.5 * m)
    end
end

# -----------------------------------------------------------------------
# nan* family
#
# Ported to Julia's own reductions rather than transliterated, EXCEPT where
# the MATLAB versions have behaviour a reasonable reimplementation would
# not reproduce. Those exceptions are called out per function.
# -----------------------------------------------------------------------

_nanslices(x::AbstractArray, dim::Integer) = eachslice(x; dims = 3 - dim)

function _nanreduce(f, x::AbstractMatrix, dim::Integer)
    dim == 1 && return reshape([f(Float64[v for v in c if !isnan(v)]) for c in eachcol(x)], 1, :)
    return reshape([f(Float64[v for v in r if !isnan(v)]) for r in eachrow(x)], :, 1)
end

_nanmean_(w) = isempty(w) ? NaN : mean(w)
_nansum_(w) = isempty(w) ? NaN : sum(w)
_nanstd_(w) = length(w) < 2 ? (isempty(w) ? NaN : 0.0) : std(w; corrected = true)
_nanvar_(w) = length(w) < 2 ? (isempty(w) ? NaN : 0.0) : var(w; corrected = true)
_nansem_(w) = isempty(w) ? NaN : _nanstd_(w) / sqrt(length(w))

"""
    nanmean(x, dim)

Mean ignoring NaN. An all-NaN slice yields `NaN`.
"""
nanmean(x::AbstractMatrix, dim::Integer = 1) = _nanreduce(_nanmean_, x, dim)

"""
    nansum(x, dim)

Sum ignoring NaN. **An all-NaN slice yields `NaN`, not zero** -- matching
`Auxiliary/nansum.m`, and differing from a plain `sum(skipmissing(...))`.
"""
nansum(x::AbstractMatrix, dim::Integer = 1) = _nanreduce(_nansum_, x, dim)

"""
    nanstd(x, dim)

Standard deviation ignoring NaN, `n-1` denominator (upstream `flag = 0`).
Upstream uses `max(count-1, 1)`, so a slice with exactly one non-NaN value
yields `0.0` rather than `NaN`. Replicated.
"""
nanstd(x::AbstractMatrix, dim::Integer = 1) = _nanreduce(_nanstd_, x, dim)

"""
    nanvar(x, dim)

Variance ignoring NaN, `n-1` denominator. Same single-observation
convention as [`nanstd`](@ref).
"""
nanvar(x::AbstractMatrix, dim::Integer = 1) = _nanreduce(_nanvar_, x, dim)

"""
    nansem(x, dim)

Standard error of the mean ignoring NaN.
"""
nansem(x::AbstractMatrix, dim::Integer = 1) = _nanreduce(_nansem_, x, dim)

"""
    nanmedian(x, dim)

Median ignoring NaN.

!!! danger "Upstream defect, replicated"
    `Auxiliary/nanmedian.m` ends with `i = find(y==0); y(i) = i + nan;`,
    which converts **any slice whose median is exactly zero into `NaN`**.
    That is a bug, not a NaN-propagation rule, but it is replicated here
    under the port's replicate-and-document policy: silently fixing it
    would make this function disagree with MATLAB on symmetric data.
    Use `median(skipmissing(...))` directly if you want the correct value.
"""
function nanmedian(x::AbstractMatrix, dim::Integer = 1)
    out = _nanreduce(w -> isempty(w) ? NaN : median(w), x, dim)
    for i in eachindex(out)
        out[i] == 0 && (out[i] = NaN)      # replicated upstream defect
    end
    return out
end

"""
    nanmax(x); nanmin(x)

`Auxiliary/nanmax.m` and `nanmin.m` are legacy shims that simply forward to
MATLAB's `max`/`min`, which already ignore NaN. Provided for API parity.
"""
nanmax(x::AbstractArray) = maximum(v for v in x if !isnan(v))
nanmin(x::AbstractArray) = minimum(v for v in x if !isnan(v))

# `nanmovavg`, `nancumsum` are NOT ported: `nancumsum` is unreferenced
# anywhere outside Auxiliary/, and `nanmovavg.m` is broken upstream -- its
# inner loop reads from the output array `y` (all NaN at that point) instead
# of the input `x`, so it returns all-NaN for every input. See PORTING_LOG.
