# =======================================================================
# lp.jl
#
# Port of v3dot1/VAR/LPoption.m, LPmakex.m, LPmodel.m -- local projections
# (Jorda 2005).
#
# This is a PARALLEL, INDEPENDENT API. It does not share the VAR/VARopt
# pipeline: LPmodel takes its own options struct and returns its own result
# type, and none of VARir/VARvd/VARhd apply to it.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

"""
    LPOptions

Options for [`LPmodel`](@ref). Immutable; update with `@set`, exactly as for
[`VAROptions`](@ref).

!!! warning "`nlag` and `const_` are absent from `LPoption.m`"
    `LPoption.m` sets only `nsteps`, `pctg`, `EXOG`, `LAGS` and `impact` —
    yet `LPmodel.m` reads `LPopt.nlag` and `LPopt.const` on its 8th and 9th
    lines. Calling `LPmodel` with an unmodified `LPoption` struct therefore
    fails in MATLAB with a bare "Reference to non-existent field 'nlag'".

    (`LPmodel.m` does contain `if ~exist('nlag','var'); error(...)` and
    `if ~exist('const','var'); const = 1; end`, but both are dead: the
    variables were already assigned from `LPopt` two lines earlier, so the
    guards can never fire.)

    This port declares both fields. `nlag` defaults to `nothing` and raises
    an informative error if left unset — preserving the *contract* that it
    is mandatory while replacing an obscure field-access failure with a
    message that says what to do. `const_` defaults to `1`, matching the
    dead guard's evident intent.

| field | default | meaning |
|---|---|---|
| `nsteps` | `40` | horizons `h = 1 … nsteps` |
| `pctg` | `95` | confidence level for the bands |
| `EXOG` | `nothing` | contemporaneous exogenous regressors |
| `LAGS` | `nothing` | variables to include as lags (drives the constant, see below) |
| `impact` | `0` | `0` one-standard-deviation shock, `1` unit shock |
| `nlag` | `nothing` | **required**; number of lags of `LAGS` |
| `const_` | `1` | deterministics passed to `VARmakexy` when building the lag block |
"""
Base.@kwdef struct LPOptions
    nsteps::Int = 40
    pctg::Float64 = 95.0
    EXOG::Union{Nothing,Matrix{Float64}} = nothing
    LAGS::Union{Nothing,Matrix{Float64}} = nothing
    impact::Int = 0
    nlag::Union{Nothing,Int} = nothing
    const_::Int = 1
end

"""
    LPoption() -> LPOptions

Default local-projection options, matching `LPoption.m` — plus the `nlag`
and `const_` fields that file omits. See [`LPOptions`](@ref).
"""
LPoption() = LPOptions()

function Base.:(==)(a::LPOptions, b::LPOptions)
    return all(f -> getfield(a, f) == getfield(b, f), fieldnames(LPOptions))
end

"""
    HorizonResult

Per-horizon regression output from [`LPmodel`](@ref), reached as
`LPout.h[hh]`.

Replaces `LPmodel.m:132`'s `eval(['aux = {''h' num2str(hh) '''};'])`, which
builds fields named `LPout.h1`, `LPout.h2`, … This is the 17th and last
`eval` site in the toolbox, and the only one that is not the per-equation
pattern handled by [`EquationResult`](@ref).

!!! warning "`bstd_HW` holds Newey-West values, not White"
    `LPmodel.m` assigns `bstd_HW` twice in consecutive statements:

    ```matlab
    LPout.(aux{1}).bstd_HW = OLSout.bstd_HW;
    LPout.(aux{1}).bstd_HW = OLSout.bstd_NW;   % overwrites the line above
    ```

    so the stored field carries the **Newey-West** standard errors and the
    White ones are lost. Replicated, with both exposed here: `bstd_HW` is
    the (mislabelled) Newey-West value MATLAB stores, and `bstd_White` is
    the value that was overwritten. Note that the `INF`/`SUP` bands are
    computed from `OLSout.bstd_HW` — the *White* errors — so the bands and
    this field do not use the same estimator.
"""
struct HorizonResult
    nobs::Int
    beta::Vector{Float64}
    tstat::Vector{Float64}
    bstd::Vector{Float64}
    bstd_HW::Vector{Float64}
    bstd_White::Vector{Float64}
    tprob::Vector{Float64}
    resid::Vector{Float64}
    yhat::Vector{Float64}
    y::Vector{Float64}
    x::Matrix{Float64}
    rsqr::Float64
    rbar::Float64
    dw::Float64
    sigma::Float64
end

"""
    LPOutput

Result of [`LPmodel`](@ref). `h` is a vector of [`HorizonResult`](@ref),
indexed by horizon.
"""
Base.@kwdef struct LPOutput
    ENDO::Matrix{Float64}
    EXOG::Union{Nothing,Matrix{Float64}}
    LAGS::Union{Nothing,Matrix{Float64}}
    S::Matrix{Float64}
    nvar::Int
    nlag::Int
    ntotcoeff::Int
    constant::Int
    h::Vector{HorizonResult}
    IR::Vector{Float64}
end

"""
    LPmakex(x, nlag)

Horizontal stack of lags `1 … nlag` of `x`, NaN-padded at the head.

Built on [`L`](@ref) with a negative argument, since in this toolbox a
*positive* argument to `L` is a lead. Matches `LPmakex.m`.
"""
function LPmakex(x::AbstractMatrix, nlag::Integer)
    nlag < 1 && return similar(x, size(x, 1), 0)
    return hcat((L(x, -ii) for ii in 1:nlag)...)
end
LPmakex(x::AbstractVector, nlag::Integer) = LPmakex(reshape(x, :, 1), nlag)

"""
    LPmodel(ENDO, S, LPopt=LPoption()) -> (IR, INF, SUP, LPout)

Local projections of `ENDO`'s **first column** on the shock series `S`,
horizon by horizon.

At horizon `h` the regression is `y[t+h-1] ~ shock[t] + EXOG[t] + lags`, so
`IR[h]` is the coefficient on the shock. `IR`, `INF` and `SUP` are
`nsteps`-vectors.

`LPopt.impact = 0` standardises `S` (z-score, `n-1` denominator) so the
response is to a one-standard-deviation shock; `impact = 1` uses `S` as is.

!!! note "Where the constant comes from"
    The per-horizon regression is run with `OLSmodel(Y, X, 0)` — **no
    constant added**. A constant enters only through the lag block, which is
    built by `VARmakexy(LAGS, nlag, const_)`. So if `LPopt.LAGS` is
    `nothing`, the projection has **no intercept at all**, whatever
    `const_` says. Upstream behaviour, preserved.

!!! warning "Bands are Normal, not bootstrapped, and use White errors"
    `INF`/`SUP` are `IR ± z * bstd_HW[1]` with `z` the Normal quantile at
    `pctg` — asymptotic, not a bootstrap, and with no Newey-West correction
    for the serial correlation that overlapping horizons induce. Since
    local-projection residuals are autocorrelated by construction for
    `h > 1`, these bands are anti-conservative at longer horizons. That is
    upstream's choice, replicated; see `LPout.h[hh].bstd_HW` for the
    Newey-West errors MATLAB computes but does not use here.

!!! warning "`EXOG`/`LAGS` row-count checks are unreachable upstream"
    `LPmodel.m` nests its `EXOG` and `LAGS` dimension checks *inside* the
    `if nobsS ~= nobsENDO` block and *after* its `error(...)` call, so they
    can never run. This port checks all three properly — a strict
    improvement that cannot change the result of any input MATLAB accepts.
"""
function LPmodel(ENDO::AbstractMatrix, S::AbstractMatrix, LPopt::LPOptions = LPoption())
    LPopt.nlag === nothing && throw(ArgumentError(
        "LPmodel: LPopt.nlag must be set. LPoption.m does not define this " *
        "field, so MATLAB fails here with a non-existent-field error; set " *
        "it explicitly, e.g. `LPopt = @set LPopt.nlag = 4`."))
    nlag = LPopt.nlag
    nlag < 1 && throw(ArgumentError(
        "Number of lags (nlag) needs to be greater or equal than 1"))

    H = LPopt.nsteps
    impact = LPopt.impact
    constant = LPopt.const_
    EXOG, LAGS = LPopt.EXOG, LPopt.LAGS

    nobsENDO = size(ENDO, 1)
    size(S, 1) == nobsENDO ||
        throw(DimensionMismatch("S has different number of observations from ENDO"))
    nvarEXOG = EXOG === nothing ? 0 : size(EXOG, 2)
    nvarLAGS = LAGS === nothing ? 0 : size(LAGS, 2)
    # These two checks are dead code upstream (see docstring); enforced here.
    nvarEXOG > 0 && size(EXOG, 1) == nobsENDO ||
        nvarEXOG == 0 ||
        throw(DimensionMismatch("EXOG has different number of observations from ENDO"))
    nvarLAGS > 0 && size(LAGS, 1) == nobsENDO ||
        nvarLAGS == 0 ||
        throw(DimensionMismatch("LAGS has different number of observations from ENDO"))

    nvarS = size(S, 2)
    ntotcoeff = constant + nvarS + nvarEXOG + nvarLAGS * nlag

    if nvarLAGS > 0
        _, auxX = VARmakexy(Matrix{Float64}(LAGS), nlag, constant)
        lags = auxX
        endo = ENDO[(1+nlag):end, 1]
        s = S[(1+nlag):end, 1]
        exog = nvarEXOG > 0 ? Matrix{Float64}(EXOG[(1+nlag):end, :]) : zeros(length(endo), 0)
    else
        lags = zeros(nobsENDO, 0)
        endo = ENDO[:, 1]
        s = S[:, 1]
        exog = nvarEXOG > 0 ? Matrix{Float64}(EXOG) : zeros(nobsENDO, 0)
    end

    if impact == 0
        shock = (s .- mean(s)) ./ std(s)          # zscore, n-1 denominator
    elseif impact == 1
        shock = s
    else
        throw(ArgumentError("Impact must be either 0 or 1"))
    end

    LHS = endo
    RHS = hcat(reshape(shock, :, 1), exog, lags)

    IR = fill(NaN, H)
    INF = fill(NaN, H)
    SUP = fill(NaN, H)
    hs = Vector{HorizonResult}(undef, H)
    conf = quantile(Normal(), 1 - (1 - LPopt.pctg / 100) / 2)

    for hh in 1:H
        Y = LHS[hh:end]
        X = RHS[1:(end-hh+1), :]
        O = OLSmodel(Y, X, 0)
        nobse = length(Y)
        tprob = tdis_prb(O.tstat, nobse - ntotcoeff)

        hs[hh] = HorizonResult(nobse, O.beta, O.tstat, O.bstd,
                               O.bstd_NW,      # upstream stores NW under bstd_HW
                               O.bstd_HW,      # the value that overwrite discarded
                               tprob, O.resid, O.yhat, Y, X,
                               O.rsqr, O.rbar, O.dw, O.sige)

        IR[hh] = O.beta[1]
        INF[hh] = IR[hh] - conf * O.bstd_HW[1]
        SUP[hh] = IR[hh] + conf * O.bstd_HW[1]
    end

    LPout = LPOutput(
        ENDO = Matrix{Float64}(ENDO),
        EXOG = EXOG === nothing ? nothing : Matrix{Float64}(EXOG),
        LAGS = LAGS === nothing ? nothing : Matrix{Float64}(LAGS),
        S = Matrix{Float64}(S),
        nvar = nvarS + nvarEXOG + nvarLAGS,
        nlag = nlag, ntotcoeff = ntotcoeff, constant = constant,
        h = hs, IR = IR,
    )
    return IR, INF, SUP, LPout
end
