# =======================================================================
# stats.jl
#
# Port of v3dot1/Stats/ -- summary statistics, moving-window transforms and
# pairwise correlations.
#
# This folder is independent of the VAR/VARopt pipeline; nothing here feeds
# VARmodel or the identification layer.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

# Every Mov* function shares the same argument validation and the same
# "force a vector to a column" reshape.
function _movcheck(DATA::AbstractArray, window::Integer)
    window <= 0 && throw(ArgumentError("Window must be positive."))
    D = DATA isa AbstractVector ? reshape(DATA, :, 1) :
        (minimum(size(DATA)) == 1 ? reshape(DATA, :, 1) : DATA)
    window > size(D, 1) &&
        throw(ArgumentError("window must not be greater than the length of DATA."))
    return D
end

# Trailing-window reductions: value at row t uses rows t-window+1 … t, with
# the first window-1 rows NaN so the result stays conformable with the input.
function _movtrailing(f, DATA::AbstractArray, window::Integer)
    D = _movcheck(DATA, window)
    nobs, nvar = size(D)
    OUT = fill(NaN, nobs, nvar)
    for row in 1:(nobs-window+1)
        OUT[row+window-1, :] = f(@view D[row:(row+window-1), :])
    end
    return OUT
end

"""
    MovAvg(DATA, window)

Trailing moving average over `window` rows. The first `window-1` rows are
`NaN`, keeping the result conformable with `DATA`.

Uses a plain mean, so a `NaN` anywhere in the window propagates. Contrast
[`MovAvgCent`](@ref), which skips NaNs.
"""
MovAvg(DATA, window::Integer) =
    _movtrailing(v -> vec(mean(v, dims = 1)), DATA, window)

"""
    MovSum(DATA, window)

Trailing moving sum. `NaN` propagates.
"""
MovSum(DATA, window::Integer) =
    _movtrailing(v -> vec(sum(v, dims = 1)), DATA, window)

"""
    MovStd(DATA, window)

Trailing moving standard deviation, `n-1` denominator. `NaN` propagates.
"""
MovStd(DATA, window::Integer) =
    _movtrailing(v -> vec(std(v, dims = 1)), DATA, window)

"""
    MovAvgCent(DATA, window)

**Centred** moving average: row `t` averages rows `t-window … t+window`.

!!! warning "`window` means something different here"
    In `MovAvg`/`MovSum`/`MovStd` the window is a *total* length; in
    `MovAvgCent` it is a *half-width*, so the average spans `2*window+1`
    rows. Upstream inconsistency, preserved.

Unlike the trailing versions this one uses [`nanmean`](@ref), so NaNs inside
the window are skipped rather than propagated — another upstream
inconsistency between otherwise-parallel functions.

Rows outside `window+1 … nobs-window` are `NaN`.
"""
function MovAvgCent(DATA, window::Integer)
    D = _movcheck(DATA, window)
    nobs, nvar = size(D)
    OUT = fill(NaN, nobs, nvar)
    for jj in (window+1):(nobs-window)
        OUT[jj, :] = vec(nanmean(Matrix(D[(jj-window):(jj+window), :]), 1))
    end
    return OUT
end

# Pearson correlation matrix, matching MATLAB's corrcoef (NaN propagates).
function _corrmat(X::AbstractMatrix)
    n = size(X, 2)
    C = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        C[i, j] = cor(view(X, :, i), view(X, :, j))
    end
    return C
end

"""
    MovCorr(X, window, i)

Trailing moving correlation of column `i` of `X` with every other column.

Returns `(nobs × nvars-1)`; row `t` uses rows `t-window+1 … t`. Rows before
`window` are `NaN`.
"""
function MovCorr(X::AbstractMatrix, window::Integer, i::Integer)
    nobs, nvars = size(X)
    OUT = fill(NaN, nobs, nvars - 1)
    idx = setdiff(1:nvars, i)
    for tt in window:nobs
        C = _corrmat(view(X, (tt-window+1):tt, :))
        OUT[tt, :] = C[i, idx]
    end
    return OUT
end

"""
    MovCorrCent(X, window, i)

Centred moving correlation of column `i` with every other column: row `t`
uses rows `t-window … t+window`. As with [`MovAvgCent`](@ref), `window` is a
half-width here.
"""
function MovCorrCent(X::AbstractMatrix, window::Integer, i::Integer)
    nobs, nvars = size(X)
    OUT = fill(NaN, nobs, nvars - 1)
    idx = setdiff(1:nvars, i)
    for tt in (window+1):(nobs-window)
        C = _corrmat(view(X, (tt-window):(tt+window), :))
        OUT[tt, :] = C[i, idx]
    end
    return OUT
end

"""
    censor(X, pct)

Set values below the `pct`-th or above the `(100-pct)`-th percentile of each
column to `NaN`.

Percentiles use MATLAB's definition (Hazen/type-5) via `_prctile`, not
Julia's default — see M1.
"""
function censor(X::AbstractMatrix, pct::Real)
    OUT = fill(NaN, size(X))
    for ii in 1:size(X, 2)
        x = Vector{Float64}(X[:, ii])
        upp = _prctile1(x, 100 - pct)
        low = _prctile1(x, pct)
        for k in eachindex(x)
            (x[k] > upp || x[k] < low) && (x[k] = NaN)
        end
        OUT[:, ii] = x
    end
    return OUT
end

"""
    CorrTable(x; pairwise=false) -> CORR

Lower-triangular correlation matrix with a unit diagonal; the upper triangle
is `NaN`.

Set `pairwise = true` to compute each pair on its own common sample rather
than dropping any row containing a `NaN`.

!!! note "Only the numeric output is ported"
    `CorrTable.m` also returns a cell-array `TABLE` and prints via `mprint`.
    The display layer is not ported — see PORTING_LOG.

!!! warning "An exactly-zero correlation is not blanked upstream"
    `CorrTable.m` blanks the upper triangle with
    `CORR(find(triu(CORR))) = NaN`, and `find` returns only **non-zero**
    entries. An upper-triangular correlation that is exactly `0.0`
    therefore survives as `0.0` instead of becoming `NaN`. Replicated.
"""
function CorrTable(x::AbstractMatrix; pairwise::Bool = false)
    n = size(x, 2)
    CORR = Matrix{Float64}(undef, n, n)
    if pairwise
        for i in 1:n, j in 1:n
            a, b = view(x, :, i), view(x, :, j)
            keep = .!(isnan.(a) .| isnan.(b))
            CORR[i, j] = any(keep) ? cor(a[keep], b[keep]) : NaN
        end
    else
        keep = .!vec(any(isnan, x, dims = 2))
        CORR = _corrmat(x[keep, :])
    end
    # Blank the upper triangle -- but only non-zero entries, as `find` does.
    for i in 1:n, j in i:n
        CORR[i, j] != 0 && (CORR[i, j] = NaN)
    end
    for i in 1:n
        CORR[i, i] = 1.0
    end
    return CORR
end

"""
    PairCorr(DATA; absolute=false) -> PC

Average pairwise correlation of each column with all the others, computed
both on levels and on log first differences.

Returns `(nvar × 2)`: column 1 levels, column 2 log differences.

!!! warning "Takes logs unconditionally"
    `PairCorr.m` computes `log(X[2:end,:]) - log(X[1:end-1,:])` with no
    check of sign, so any non-positive entry yields `NaN` (MATLAB) — a
    complex number would be produced for negatives were it not for the
    subsequent correlation. Only pass strictly positive levels. Replicated.
"""
function PairCorr(DATA::AbstractMatrix; absolute::Bool = false)
    X = Matrix{Float64}(DATA)
    # Upstream blanks any column whose FIRST observation is NaN.
    for j in 1:size(X, 2)
        isnan(X[1, j]) && (X[:, j] .= NaN)
    end
    d = log.(X[2:end, :]) .- log.(X[1:(end-1), :])

    X_corr = _corrmat(X)
    x_corr = _corrmat(d)
    COR1 = absolute ? abs.(X_corr) : X_corr
    COR2 = absolute ? abs.(x_corr) : x_corr

    nans = isnan.(X_corr)
    n = Float64.(vec(sum(.!nans, dims = 1)))
    n[n.==0] .= NaN

    lvl = (vec(nansum(COR1, 1)) .- 1) ./ (n .- 1)
    dif = (vec(nansum(COR2, 1)) .- 1) ./ (n .- 1)
    return hcat(lvl, dif)
end

"""
    CrossCorrelation(X, Y, lag) -> (OUT, lags)

Cross-correlations of the `T×1` series `X` with each column of `Y` at leads
and lags `-lag … +lag`.

`OUT` is `(2*lag+1) × size(Y,2)`; `lags` is the corresponding index vector.
The plotting branch of `CrossCorrelation.m` is not ported.
"""
function CrossCorrelation(X::AbstractVecOrMat, Y::AbstractMatrix, lag::Integer)
    x = vec(X)
    a = length(x)
    size(Y, 1) == a ||
        throw(DimensionMismatch("vectors must have the same length"))
    lag >= a - 2 &&
        throw(ArgumentError("too few observations. Reduce the number of lags"))
    d = size(Y, 2)

    OUT = Matrix{Float64}(undef, 2 * lag + 1, d)
    m = 1
    for jj in 1:lag
        for ii in 1:d
            OUT[m, ii] = cor(x[1:(end-lag-1+jj)], Y[(1+lag+1-jj):end, ii])
        end
        m += 1
    end
    for ii in 1:d
        OUT[m, ii] = cor(x, Y[:, ii])
    end
    m += 1
    for jj in 1:lag
        for ii in 1:d
            OUT[m, ii] = cor(x[(1+jj):end], Y[1:(end-jj), ii])
        end
        m += 1
    end
    return OUT, collect(-lag:lag)
end
