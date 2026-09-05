# =======================================================================
# utils.jl
#
# Port of v3dot1/Utils/ plus the lag/design-matrix helpers that live in
# v3dot1/VAR/ (VARmakexy, VARmakelags, L).
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

# -----------------------------------------------------------------------
# Trivial MATLAB shims -- internal only.
#
# `cols`, `rows`, `columns` and `vec` exist upstream only because MATLAB
# lacks concise equivalents. Julia has `size(X,2)`, `size(X,1)` and
# `Base.vec`. They are kept as unexported internals so that call sites can
# stay visually close to the MATLAB source without polluting the namespace
# of a package whose users will reasonably want `rows` and `cols` as their
# own variable names.
#
# NB: upstream `Utils/columns.m` declares `function c = cols(DATA)` -- the
# function name inside the file does not match the filename. MATLAB
# dispatches on the filename so `columns(x)` works anyway.
# -----------------------------------------------------------------------
_cols(X::AbstractArray) = size(X, 2)
_rows(X::AbstractArray) = size(X, 1)

function _vec(X::AbstractArray, dim::Integer = 1)
    dim == 1 && return vec(X)
    dim == 2 && return reshape(X, 1, :)
    throw(ArgumentError("input \"dim\" has wrong value"))
end

"""
    Num2NaN(DATA)

Replace the sentinel `123456789` with `NaN`. Matches `Utils/Num2NaN.m`.
"""
Num2NaN(DATA::AbstractArray) = map(v -> v == 123456789 ? convert(float(typeof(v)), NaN) : float(v), DATA)

"""
    NaN2Num(DATA)

Replace `NaN` with the sentinel `123456789`. Matches `Utils/NaN2Num.m`.
"""
NaN2Num(DATA::AbstractArray) = map(v -> isnan(v) ? 123456789.0 : float(v), DATA)

"""
    CommonSample(DATA, dim=1) -> (OUT, fo, lo)

Drop every row (`dim=1`) or column (`dim=2`) of `DATA` containing a `NaN`.

Returns the trimmed data together with two counts. The three-value return
is load-bearing: `fo`/`lo` are how callers trim an external instrument to
the same sample the VAR was estimated on.

`fo` is the number of leading all-NaN rows (columns) and `lo` the number of
trailing ones. An empty input is returned unchanged with `fo = lo = 0`.

!!! note "Upstream v4.0 fixes adopted"
    In VAR Toolbox 3.1 the loop computing `lo` never stopped at the first
    non-NaN from the end, so `lo` counted **every** NaN row outside the
    leading run rather than the trailing run; and an empty input raised an
    index error. Both were fixed in v4.0 (`Utils/CommonSample.m`:
    `else; break` in the trailing loop; early return on empty input). This
    port follows v4.0, so `lo` is now the documented trailing count.
"""
function CommonSample(DATA::AbstractMatrix, dim::Integer = 1)
    fo = 0
    lo = 0
    isempty(DATA) && return DATA, fo, lo             # v4.0: empty in, empty out
    if dim == 1
        temp = vec(sum(DATA, dims = 2))
        n = length(temp)
        if n >= 1 && isnan(temp[1])
            ii = 1
            while isnan(temp[ii])
                fo += 1
                ii += 1
                ii > n && break
            end
        end
        for ii in 1:(size(DATA, 1)-fo)
            isnan(temp[end+1-ii]) || break           # v4.0: stop at first non-NaN
            lo += 1
        end
        keep = .!vec(any(isnan, DATA, dims = 2))
        OUT = DATA[keep, :]
    else
        temp = vec(sum(DATA, dims = 1))
        n = length(temp)
        if n >= 1 && isnan(temp[1])
            ii = 1
            while isnan(temp[ii])
                fo += 1
                ii += 1
                ii > n && break
            end
        end
        for ii in 1:(size(DATA, 2)-fo)
            isnan(temp[end+1-ii]) || break           # v4.0: stop at first non-NaN
            lo += 1
        end
        keep = .!vec(any(isnan, DATA, dims = 1))
        OUT = DATA[:, keep]
    end
    return OUT, fo, lo
end

"""
    getqr(a)

QR factorisation returning an orthonormal `Q` with the sign of each column
fixed so that `diag(R) >= 0`. Matches `Utils/getqr.m`.

The sign fix is what makes the sign-restriction rotations reproducible; do
not replace with a bare `qr(a).Q`.
"""
function getqr(a::AbstractMatrix)
    F = qr(a)
    q = Matrix(F.Q)
    r = F.R
    for i in 1:size(q, 1)
        if i <= size(r, 1) && i <= size(r, 2) && r[i, i] < 0
            @views q[:, i] .= .-q[:, i]
        end
    end
    return q
end

"""
    OrthNorm(n; rng=Random.default_rng(), draws=nothing)

A random orthonormal `n x n` matrix: `getqr(randn(n,n))`. Matches
`VAR/OrthNorm.m`.

Pass `draws` (an `n x n` matrix) to supply the underlying normal draws
explicitly -- this is how MATLAB-generated draws are injected for
golden-master testing. See M2 in the porting brief.
"""
function OrthNorm(n::Integer; rng::AbstractRNG = Random.default_rng(), draws = nothing)
    a = draws === nothing ? randn(rng, n, n) : draws
    return getqr(a)
end

"""
    trimr(x, n1, n2)

Drop the first `n1` and last `n2` rows. Matches `Auxiliary/trimr.m`,
including its error when the trim exceeds the data.
"""
function trimr(x::AbstractMatrix, n1::Integer, n2::Integer)
    n = size(x, 1)
    (n1 + n2) >= n && throw(ArgumentError("Attempting to trim too much in trimr"))
    return x[(n1+1):(n-n2), :]
end
trimr(x::AbstractVector, n1::Integer, n2::Integer) = vec(trimr(reshape(x, :, 1), n1, n2))

"""
    lag(x, n=1, v=0.0)

Shift `x` down by `n` rows, padding the head with `v`. Matches
`Auxiliary/lag.m`, including its `n < 1 -> empty` behaviour.

Unexported: `lag` is too generic a name to claim in a user's namespace.
Reach it as `VARToolbox.lag`.
"""
function lag(x::AbstractMatrix, n::Integer = 1, v::Real = 0.0)
    n < 1 && return similar(x, 0, size(x, 2))
    zt = fill(float(v), n, _cols(x))
    return vcat(zt, trimr(x, 0, n))
end

"""
    L(X, nlag)

Lead/lag operator, padding with `NaN`. Matches `VAR/L.m`.

!!! warning "The sign convention is inverted from what the name suggests"
    A **positive** `nlag` produces a **lead**; a negative `nlag` produces a
    lag. This is upstream's convention (`if nlag > 0 % this is lead`),
    preserved deliberately.
"""
function L(X::AbstractMatrix, nlag::Integer)
    if nlag > 0                                   # lead
        zt = fill(NaN, nlag, _cols(X))
        return vcat(trimr(X, nlag, 0), zt)
    elseif nlag < 0                               # lag
        zt = fill(NaN, abs(nlag), _cols(X))
        return vcat(zt, trimr(X, 0, abs(nlag)))
    else
        return copy(X)
    end
end

"""
    seqa(a, b, c)

Column vector of `c` points starting at `a` in steps of `b`. Matches
`Auxiliary/seqa.m`.
"""
seqa(a::Real, b::Real, c::Integer) = collect(range(a; step = b, length = c))

"""
    ptrend(p, nobs)

Constant plus polynomial trends of order `1:p`, with time normalised by
`nobs` (`t = (1:nobs)/nobs`). Matches `Auxiliary/ptrend.m`.
"""
function ptrend(p::Integer, nobs::Integer)
    u = ones(nobs, 1)
    p <= 0 && return u
    t = collect(1:nobs) ./ nobs
    timep = hcat((t .^ m for m in 1:p)...)
    return hcat(u, timep)
end

# -----------------------------------------------------------------------
# VAR design matrices
# -----------------------------------------------------------------------

"""
    VARmakelags(DATA, lag) -> OUT

Stack `DATA` with its first `lag` lags: `[y_t, y_{t-1}, ..., y_{t-lag}]`,
losing the first `lag` rows. Matches `VAR/VARmakelags.m`.
"""
function VARmakelags(DATA::AbstractMatrix, lag::Integer)
    nobs = size(DATA, 1)
    blocks = [DATA[(jj+1):(nobs-lag+jj), :] for jj in 0:(lag-1)]
    OUT = isempty(blocks) ? similar(DATA, nobs, 0) : hcat(reverse(blocks)...)
    return hcat(DATA[(lag+1):end, :], OUT)
end

"""
    VARmakexy(DATA, lags, const_) -> (Y, X)

Build the VAR regressand and regressor matrices. Matches `VAR/VARmakexy.m`.

`const_` selects the deterministics:

| value | regressors prepended |
|---|---|
| 0 | none |
| 1 | constant |
| 2 | constant, linear trend |
| 3 | constant, linear trend, squared trend |

The trend is the raw `1:T` sequence, **not** normalised by `T` (unlike
[`ptrend`](@ref)).

Named `const_` because `const` is reserved in Julia.
"""
function VARmakexy(DATA::AbstractMatrix, lags::Integer, const_::Integer)
    nobs = size(DATA, 1)
    Y = DATA[(lags+1):end, :]
    blocks = [DATA[(jj+1):(nobs-lags+jj), :] for jj in 0:(lags-1)]
    Xlag = isempty(blocks) ? similar(DATA, nobs - lags, 0) : hcat(reverse(blocks)...)

    T = nobs - lags
    if const_ == 0
        X = Xlag
    elseif const_ == 1
        X = hcat(ones(T), Xlag)
    elseif const_ == 2
        trend = collect(1.0:T)
        X = hcat(ones(T), trend, Xlag)
    elseif const_ == 3
        trend = collect(1.0:T)
        X = hcat(ones(T), trend, trend .^ 2, Xlag)
    else
        throw(ArgumentError("VARmakexy: const must be 0, 1, 2 or 3 (got $const_)"))
    end
    return Y, X
end
