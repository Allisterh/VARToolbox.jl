# =======================================================================
# identification.jl
#
# Port of v3dot1/VAR/VARir.m, VARvd.m, VARhd.m.
#
# Identification is LATE-BOUND: VARmodel performs none of it. The scheme is
# selected by VARopt.ident and applied here. Each function returns an
# updated model as its second output, with B / PSI / Fp filled in.
#
# VAR Toolbox 3.1 -> VARToolbox.jl
# Original MATLAB: Ambrogio Cesa-Bianchi
# =======================================================================

# -----------------------------------------------------------------------
# Wold multipliers, shared by VARir and VARvd
# -----------------------------------------------------------------------

# Fp[:,:,i] is the i-th lag coefficient block; PSI[:,:,i] the i-th Wold
# multiplier. Matches VARir.m:64-79 and VARvd.m (identical code).
#
# NB on the size of Fp. MATLAB preallocates `zeros(nvar,nvar,nlag)` but its
# else-branch then assigns `Fp(:,:,ii) = zeros(nvar,nvar)` for ii > nlag,
# which *grows* the array to nsteps slices. The Wold recursion below indexes
# Fp[:,:,jj] for jj up to nsteps-1, so those zero slices are load-bearing --
# allocating only nlag of them runs off the end. Hence max(nlag, nsteps).
function _wold(VAR::VARModel, nsteps::Integer)
    nvar, nlag = VAR.nvar, VAR.nlag
    Fp = zeros(nvar, nvar, max(nlag, nsteps))
    I0 = VAR.constant + 1
    for ii in 1:nsteps
        if ii <= nlag
            Fp[:, :, ii] = VAR.F[:, I0:(I0+nvar-1)]
        end
        I0 += nvar
    end
    PSI = zeros(nvar, nvar, nsteps)
    PSI[:, :, 1] = Matrix{Float64}(I, nvar, nvar)
    for ii in 2:nsteps
        aux = zeros(nvar, nvar)
        for jj in 1:(ii-1)
            aux += PSI[:, :, ii-jj] * Fp[:, :, jj]
        end
        PSI[:, :, ii] = aux
    end
    return PSI, Fp
end

const _IDENT_HELP = """
Choose one of the following options:
- short: zero contemporaneous restrictions
- long:  zero long-run restrictions
- sign:  sign restrictions
- iv:    external instrument
- exog:  exogenous variable"""

# The `short` / `long` / `sign` branches, which all three of VARir, VARvd and
# VARhd share verbatim.
function _B_common(VAR::VARModel, ident::AbstractString)
    if ident == "short"
        C = cholesky(Symmetric(VAR.sigma); check = false)
        issuccess(C) || error("VCV is not positive definite")
        return Matrix(C.U)'
    elseif ident == "long"
        Finf_big = inv(Matrix{Float64}(I, size(VAR.Fcomp)...) - VAR.Fcomp)
        Finf = Finf_big[1:VAR.nvar, 1:VAR.nvar]
        D = Matrix(cholesky(Symmetric(Finf * VAR.sigma * transpose(Finf))).U)'
        return Finf \ D
    elseif ident == "sign"
        isempty(VAR.B) && error(
            "You need to provide the B matrix with SR.m and/or SignRestrictions.m")
        return VAR.B
    end
    return nothing
end

# -----------------------------------------------------------------------
# VARir
# -----------------------------------------------------------------------

"""
    VARir(VAR, VARopt) -> (IR, VAR)

Impulse responses.

`IR` is `(nsteps × nvar × nshocks)` — horizon, then **response variable**,
then **shock**.

!!! warning "`VD` uses the opposite convention"
    [`VARvd`](@ref) returns `(nsteps × nshocks × nvar)` — shock and variable
    swapped relative to `IR`. This is upstream behaviour, verified
    numerically (a variance decomposition sums to 100 across `VD`'s
    *second* index). Do not assume the two arrays are indexed alike.

Returns an updated model as its second output, carrying `B`, `PSI` and `Fp`.
Capture it — the input model is unchanged:

```julia
IR, VAR = VARir(VAR, VARopt)
```

# Identification

Set `VARopt.ident` before calling. `VARir` supports all five schemes:

| `ident` | meaning |
|---|---|
| `"short"` | Cholesky, zero contemporaneous restrictions |
| `"long"` | zero long-run restrictions |
| `"sign"` | uses a `B` already placed by [`SR`](@ref) |
| `"iv"` | external instrument |
| `"exog"` | first exogenous variable's coefficient column |

!!! warning "`ident = \"iv\"` instruments the FIRST column of ENDO"
    Variable order is load-bearing. The variable to be instrumented must be
    ordered first in `ENDO`; `VAR.IV` must be set. Reordering `ENDO` changes
    the result.

!!! note "`ident = \"sign\"` requires `B` to be set already"
    Setting `ident = "sign"` by hand is not a working path — it raises.
    Go through [`SR`](@ref), which places `B` for you.
"""
function VARir(VAR::VARModel, VARopt::VAROptions)
    ident = VARopt.ident
    if ident == "iv" && isempty(VAR.IV)
        error("You need to provide the data for the instrument in VAR (IV)")
    end

    nsteps, impact, shut, recurs = VARopt.nsteps, VARopt.impact, VARopt.shut, VARopt.recurs
    nvar, nlag = VAR.nvar, VAR.nlag
    Fcomp = copy(VAR.Fcomp)
    IR = fill(NaN, nsteps, nvar, nvar)

    PSI, Fp = _wold(VAR, nsteps)

    FirstStage = nothing
    sigma_b = zeros(0, 0)
    Biv = zeros(0, 0)

    B = _B_common(VAR, ident)
    if B === nothing
        if ident == "iv"
            up = VAR.resid[:, 1]
            uq = VAR.resid[:, 2:end]

            aux, _, _ = CommonSample(hcat(up, VAR.IV[(nlag+1):end, :]))
            p = aux[:, 1]
            q = uq[(end-length(p)+1):end, :]
            pq = hcat(p, q)
            Z = aux[:, 2:end]

            FirstStage = OLSmodel(p, Z)
            p_hat = FirstStage.yhat

            Bivv = zeros(nvar)
            Bivv[1] = 1.0                        # impact response normalised to 1
            sqsp = zeros(size(q, 2))
            for ii in 2:nvar
                SecondStage = OLSmodel(q[:, ii-1], reshape(p_hat, :, 1))
                Bivv[ii] = SecondStage.beta[2]
                sqsp[ii-1] = SecondStage.beta[2]
            end

            # Rescale to the true shock size -- Gertler & Karadi (2015) fn 4.
            dm = pq .- transpose(mean(pq, dims = 1))'
            sigma_b = (1 / (size(pq, 1) - VAR.ntotcoeff)) .* (transpose(dm) * dm)
            s21s11 = sqsp
            S11 = sigma_b[1, 1]
            S21 = sigma_b[2:end, 1]
            S22 = sigma_b[2:end, 2:end]
            Q = s21s11 * S11 * transpose(s21s11) -
                (S21 * transpose(s21s11) + s21s11 * transpose(S21)) + S22
            r = S21 - s21s11 .* S11
            sp = sqrt(S11 - (transpose(r)*(Q \ r))[1])
            Bivv .*= sp
            Biv = reshape(Bivv, :, 1)
            B = zeros(nvar, nvar)
            B[:, 1] = Bivv
        elseif ident == "exog"
            B = zeros(nvar, nvar)
            B[:, 1] = VAR.F[:, VAR.constant+VAR.ncoeff+1]
        else
            error("Identification incorrectly specified.\n" * _IDENT_HELP)
        end
    end

    # ---- responses ------------------------------------------------------
    for mm in 1:nvar
        shut != 0 && (Fcomp[shut, :] .= 0)
        response = zeros(nvar, nsteps)
        impulse = zeros(nvar)
        if impact == 0
            impulse[mm] = 1.0                    # one standard deviation
        elseif impact == 1
            impulse[mm] = 1 / B[mm, mm]          # unit shock
        else
            throw(ArgumentError("Impact must be either 0 or 1"))
        end
        response[:, 1] = B * impulse
        shut != 0 && (response[shut, 1] = 0)
        if recurs == "wold"
            for kk in 2:nsteps
                response[:, kk] = PSI[:, :, kk] * B * impulse
            end
        elseif recurs == "comp"
            for kk in 2:nsteps
                FcompN = Fcomp^(kk - 1)
                response[:, kk] = FcompN[1:nvar, 1:nvar] * B * impulse
            end
        end
        IR[:, :, mm] = transpose(response)
    end

    VAR2 = @set VAR.PSI = PSI
    VAR2 = @set VAR2.Fp = Fp
    VAR2 = @set VAR2.B = B
    if ident == "iv"
        VAR2 = @set VAR2.FirstStage = FirstStage
        VAR2 = @set VAR2.sigma_b = sigma_b
        VAR2 = @set VAR2.Biv = Biv
    end
    return IR, VAR2
end

# -----------------------------------------------------------------------
# VARvd
# -----------------------------------------------------------------------

"""
    VARvd(VAR, VARopt) -> (VD, VAR)

Forecast error variance decomposition, in percent.

`VD` is `(nsteps × nshocks × nvar)`.

!!! warning "This is NOT the same layout as `IR`"
    [`VARir`](@ref) returns `(nsteps × nvar × nshocks)`; `VD` swaps the last
    two. Verified numerically: `sum(VD[h, :, k])` is 100 for every horizon
    `h` and variable `k`, i.e. contributions sum across the **second**
    index. Upstream inconsistency, preserved.

!!! warning "Only `short`, `long` and `sign` are supported"
    `VARvd` raises on `ident = "iv"` (upstream message: "not available yet")
    and on `ident = "exog"` — the `exog` branch exists only in
    [`VARir`](@ref), so `exog` falls through to the generic
    "Identification incorrectly specified" error. Both are upstream
    behaviour, not port limitations.
"""
function VARvd(VAR::VARModel, VARopt::VAROptions)
    ident = VARopt.ident
    ident == "iv" && error(
        "Forecast error variance decomposition not available with " *
        "external instruments identification (iv)")

    nsteps = VARopt.nsteps
    nvar = VAR.nvar
    sigma = VAR.sigma
    VD = zeros(nsteps, nvar, nvar)
    SE = zeros(nsteps, nvar)

    PSI, Fp = _wold(VAR, nsteps)

    B = _B_common(VAR, ident)
    B === nothing && error("Identification incorrectly specified.\n" * _IDENT_HELP)

    MSE = zeros(nvar, nvar, nsteps)
    MSE_shock = zeros(nvar, nvar, nsteps)
    for ii in 1:nvar                              # loop over shocks
        MSE[:, :, 1] = sigma
        for nn in 2:nsteps
            MSE[:, :, nn] = MSE[:, :, nn-1] + PSI[:, :, nn] * sigma * transpose(PSI[:, :, nn])
        end
        MSE_shock[:, :, 1] = B[:, ii] * transpose(B[:, ii])
        for nn in 2:nsteps
            MSE_shock[:, :, nn] = MSE_shock[:, :, nn-1] +
                                  PSI[:, :, nn] * MSE_shock[:, :, 1] * transpose(PSI[:, :, nn])
        end
        FECD = MSE_shock ./ MSE
        for nn in 1:nsteps
            for kk in 1:nvar
                VD[nn, ii, kk] = 100 * FECD[kk, kk, nn]
            end
            SE[nn, :] = sqrt.(diag(MSE[:, :, nn]))
        end
    end

    VAR2 = @set VAR.PSI = PSI
    VAR2 = @set VAR2.Fp = Fp
    VAR2 = @set VAR2.B = B
    return VD, VAR2
end

# -----------------------------------------------------------------------
# VARhd
# -----------------------------------------------------------------------

"""
    HDOutput

Historical decomposition components from [`VARhd`](@ref). Every component is
NaN-padded at the front so that all of them align with the original sample.

| field | shape |
|---|---|
| `shock` | `(nobs+nlag) × nshocks × nvar` |
| `init`, `constant`, `trend`, `trend2`, `endo` | `(nobs+nlag) × nvar` |
| `exo` | `(nobs+nlag) × nvar × nvar_ex` |

`constant` is MATLAB's `HD.const` (reserved word in Julia).

!!! note "`init` pads with `nlag-1` NaNs, the others with `nlag`"
    `init` keeps its first column and pads with `nlag-1`; every other
    component drops its first column and pads with `nlag`. Both end up with
    `nobs+nlag` rows, so they align — but the leading NaN counts differ by
    one. Upstream, preserved.
"""
Base.@kwdef struct HDOutput
    shock::Array{Float64,3}
    init::Matrix{Float64}
    constant::Matrix{Float64}
    trend::Matrix{Float64}
    trend2::Matrix{Float64}
    exo::Array{Float64,3}
    endo::Matrix{Float64}
end

"""
    VARhd(VAR, VARopt) -> (HD, VAR)

Historical decomposition: additive contributions of each structural shock,
the initial condition, the deterministics, and the exogenous variables.

`HD.shock` is `(nobs+nlag) × nshocks × nvar` — the same shock-then-variable
ordering as [`VARvd`](@ref), not [`VARir`](@ref)'s.

!!! warning "Only `short`, `long` and `sign` are supported"
    Raises on `ident = "iv"` and on `ident = "exog"`, matching upstream.

!!! warning "Exogenous contributions ignore lags of the exogenous variables"
    `VARhd.m` builds the exogenous contribution from
    `F[:, nvar*nlag+const+ii]` only — the *contemporaneous* coefficient on
    exogenous variable `ii`. When `nlag_ex > 0` the lagged exogenous
    coefficients are silently omitted, so `HD.exo` does not account for
    their contribution and `HD.endo` will not reconstruct the data exactly.
    Upstream defect, replicated.
"""
function VARhd(VAR::VARModel, VARopt::VAROptions)
    ident = VARopt.ident
    ident == "iv" && error(
        "Historical decomposition not available with " *
        "external instruments identification (iv)")

    B = _B_common(VAR, ident)
    B === nothing && error("Identification incorrectly specified.\n" * _IDENT_HELP)

    constant = VAR.constant
    F = Matrix(transpose(VAR.Ft))
    nvar, nvar_ex = VAR.nvar, VAR.nvar_ex
    nlag = VAR.nlag
    nvarXeq = nvar * nlag
    Fcomp = VAR.Fcomp
    Y = VAR.Y
    X = VAR.X[:, (1+constant):(nvarXeq+constant)]
    nobs = size(Y, 1)

    eps_ = B \ transpose(VAR.resid)               # structural errors

    B_big = zeros(nvarXeq, nvar)
    B_big[1:nvar, :] = B
    Icomp = hcat(Matrix{Float64}(I, nvar, nvar), zeros(nvar, (nlag - 1) * nvar))

    # ---- contribution of each shock ------------------------------------
    HDshock_big = zeros(nlag * nvar, nobs + 1, nvar)
    HDshock = zeros(nvar, nobs + 1, nvar)
    for j in 1:nvar
        eps_big = zeros(nvar, nobs + 1)
        eps_big[j, 2:end] = eps_[j, :]
        for i in 2:(nobs+1)
            HDshock_big[:, i, j] = B_big * eps_big[:, i] + Fcomp * HDshock_big[:, i-1, j]
            HDshock[:, i, j] = Icomp * HDshock_big[:, i, j]
        end
    end

    # ---- initial condition ---------------------------------------------
    HDinit_big = zeros(nlag * nvar, nobs + 1)
    HDinit = zeros(nvar, nobs + 1)
    HDinit_big[:, 1] = X[1, :]
    HDinit[:, 1] = Icomp * HDinit_big[:, 1]
    for i in 2:(nobs+1)
        HDinit_big[:, i] = Fcomp * HDinit_big[:, i-1]
        HDinit[:, i] = Icomp * HDinit_big[:, i]
    end

    # ---- deterministics -------------------------------------------------
    HDconst_big = zeros(nlag * nvar, nobs + 1)
    HDconst = zeros(nvar, nobs + 1)
    if constant > 0
        CC = zeros(nlag * nvar)
        CC[1:nvar] = F[:, 1]
        for i in 2:(nobs+1)
            HDconst_big[:, i] = CC + Fcomp * HDconst_big[:, i-1]
            HDconst[:, i] = Icomp * HDconst_big[:, i]
        end
    end

    HDtrend_big = zeros(nlag * nvar, nobs + 1)
    HDtrend = zeros(nvar, nobs + 1)
    if constant > 1
        TT = zeros(nlag * nvar)
        TT[1:nvar] = F[:, 2]
        for i in 2:(nobs+1)
            HDtrend_big[:, i] = TT .* (i - 1) + Fcomp * HDtrend_big[:, i-1]
            HDtrend[:, i] = Icomp * HDtrend_big[:, i]
        end
    end

    HDtrend2_big = zeros(nlag * nvar, nobs + 1)
    HDtrend2 = zeros(nvar, nobs + 1)
    if constant > 2
        TT2 = zeros(nlag * nvar)
        TT2[1:nvar] = F[:, 3]
        for i in 2:(nobs+1)
            HDtrend2_big[:, i] = TT2 .* ((i - 1)^2) + Fcomp * HDtrend2_big[:, i-1]
            HDtrend2[:, i] = Icomp * HDtrend2_big[:, i]
        end
    end

    # ---- exogenous ------------------------------------------------------
    HDexo = zeros(nvar, nobs + 1, max(nvar_ex, 1))
    if nvar_ex > 0
        for ii in 1:nvar_ex
            HDexo_big = zeros(nlag * nvar, nobs + 1)
            VARexo = VAR.X_EX[:, ii]
            EXOcol = zeros(nlag * nvar)
            EXOcol[1:nvar] = F[:, nvar*nlag+constant+ii]
            for i in 2:(nobs+1)
                HDexo_big[:, i] = EXOcol .* VARexo[i-1] + Fcomp * HDexo_big[:, i-1]
                HDexo[:, i, ii] = Icomp * HDexo_big[:, i]
            end
        end
    end

    HDendo = HDinit + HDconst + HDtrend + HDtrend2 +
             dropdims(sum(HDexo, dims = 3), dims = 3) +
             dropdims(sum(HDshock, dims = 3), dims = 3)

    # ---- assemble, NaN-padded to the original sample --------------------
    shock = zeros(nobs + nlag, nvar, nvar)
    for i in 1:nvar, j in 1:nvar
        shock[:, j, i] = vcat(fill(NaN, nlag), HDshock[i, 2:end, j])
    end
    pad(M) = vcat(fill(NaN, nlag, nvar), transpose(M[:, 2:end]))

    exo = zeros(nobs + nlag, nvar, nvar_ex)
    for i in 1:nvar_ex
        exo[:, :, i] = pad(HDexo[:, :, i])
    end

    HD = HDOutput(
        shock = shock,
        init = vcat(fill(NaN, nlag - 1, nvar), transpose(HDinit)),
        constant = pad(HDconst),
        trend = pad(HDtrend),
        trend2 = pad(HDtrend2),
        exo = exo,
        endo = pad(HDendo),
    )

    VAR2 = @set VAR.B = B
    return HD, VAR2
end
