# =======================================================================
# replications.jl
#
# The four published replications shipped with the MATLAB toolbox,
# rewritten against VARToolbox.jl:
#
#   Stock & Watson  (2001, JEP)   -- recursive (Cholesky) identification
#   Blanchard & Quah (1989, AER)  -- long-run restrictions
#   Gertler & Karadi (2015, AEJ)  -- external instrument
#   Uhlig           (2005, JME)   -- sign restrictions
#
# These read the SAME .xlsx fixtures as v3dot1/Replic/*/GO_*.m, so there is
# one data source for both languages. Each returns a NamedTuple rather than
# plotting, which lets test/test_phase6.jl assert against MATLAB output
# instead of comparing figures by eye.
#
# Requires XLSX.jl (a test/example dependency, not a package dependency --
# VARToolbox itself does no file I/O).
#
# Usage:
#     include("replications.jl")
#     sw = stock_watson_2001()
#     sw.IR    # (nsteps x nvar x nshocks)
# =======================================================================

using VARToolbox
using XLSX
using Random
using LinearAlgebra

const REPLIC_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "v3dot1", "Replic"))

"""
    load_replication_data(relpath) -> (data, vnames, vnames_long, dates)

Read one of the toolbox's `.xlsx` fixtures.

The sheets are laid out with long names on row 1, short names on row 2, and
data from row 3 with dates in column 1 — the layout `xlsread` + the
`xlstext(2,2:end)` idiom in the `GO_*.m` scripts assumes.

Missing values are stored as the sentinel `123456789` and converted to `NaN`
here, matching `Num2NaN`. Skipping that step is not cosmetic: the
Gertler-Karadi instrument is only populated over part of the sample, and
leaving the sentinels in makes its first-stage regression numerically
singular.
"""
function load_replication_data(relpath::AbstractString)
    path = joinpath(REPLIC_DIR, relpath)
    sheet = XLSX.readxlsx(path)["Sheet1"]
    raw = sheet[:]
    vnames_long = String[string(x) for x in raw[1, 2:end]]
    vnames = String[string(x) for x in raw[2, 2:end]]
    dates = String[string(x) for x in raw[3:end, 1]]
    body = raw[3:end, 2:end]
    data = Matrix{Float64}(undef, size(body)...)
    for i in eachindex(body)
        v = body[i]
        data[i] = (v === missing) ? NaN : Float64(v)
    end
    data = Num2NaN(data)          # 123456789 sentinel -> NaN
    return data, vnames, vnames_long, dates
end

col(data, vnames, name) = data[:, findfirst(==(name), vnames)]

# -----------------------------------------------------------------------

"""
    stock_watson_2001(; ndraws=1000, band_src=nothing, vd_src=nothing)

Stock & Watson (2001, JEP), Figure 1 and Table 1.B: a trivariate VAR in
unemployment, inflation and the federal funds rate, 4 lags, recursive
identification, 24-step horizon.

Returns `(IR, B, VD, bands..., table1B)`. `table1B` reproduces the paper's
variance-decomposition table at horizons 1, 4, 8 and 12.
"""
function stock_watson_2001(; ndraws::Int = 1000,
                           band_src = nothing, vd_src = nothing,
                           rng = Random.default_rng())
    X, _, _, _ = load_replication_data(joinpath("SW2001", "SW2001_Data.xlsx"))

    VAR, VARopt = VARmodel(X, 4, 1)
    VARopt = @set VARopt.nsteps = 24
    VARopt = @set VARopt.ident = "short"
    VARopt = @set VARopt.ndraws = ndraws

    IR, VAR = VARir(VAR, VARopt)
    IRinf, IRsup, IRmed, IRbar = VARirband(VAR, VARopt; rng = rng, src = band_src)

    VD, VAR = VARvd(VAR, VARopt)
    VDinf, VDsup, VDmed, VDbar = VARvdband(VAR, VARopt; rng = rng, src = vd_src)

    # Table 1.B: rows are (variable, horizon), columns are shocks.
    # NB: VD is indexed [step, shock, variable] -- not IR's layout.
    table1B = vcat(
        (VD[h, :, 1]' for h in (1, 4, 8, 12))...,
        (VD[h, :, 2]' for h in (1, 4, 8, 12))...,
        (VD[h, :, 3]' for h in (1, 4, 8, 12))...,
    )

    return (; IR, B = VAR.B, VD, IRinf, IRsup, IRmed, IRbar,
            VDinf, VDsup, VDmed, VDbar, table1B, VAR, VARopt)
end

"""
    blanchard_quah_1989(; ndraws=1000, band_src=nothing)

Blanchard & Quah (1989, AER), Figure 2: a bivariate VAR in GDP growth and
unemployment, 8 lags, **long-run** restrictions, 40-step horizon.

The paper's figure plots the *cumulated* GDP response, and flips the sign of
the demand shock; both are reproduced in the returned `fig2_*` series.
"""
function blanchard_quah_1989(; ndraws::Int = 1000, band_src = nothing,
                             rng = Random.default_rng())
    X, _, _, _ = load_replication_data(joinpath("BQ1989", "BQ1989_Data.xlsx"))

    VAR, VARopt = VARmodel(X, 8, 1)
    VARopt = @set VARopt.nsteps = 40
    VARopt = @set VARopt.ident = "long"
    VARopt = @set VARopt.ndraws = ndraws

    IR, VAR = VARir(VAR, VARopt)
    IRinf, IRsup, IRmed, IRbar = VARirband(VAR, VARopt; rng = rng, src = band_src)

    return (; IR, B = VAR.B, IRinf, IRsup, IRmed, IRbar,
            fig2_supply_gdp = cumsum(IR[:, 1, 1]),
            fig2_supply_u = IR[:, 2, 1],
            fig2_demand_gdp = cumsum(-IR[:, 1, 2]),
            fig2_demand_u = -IR[:, 2, 2],
            VAR, VARopt)
end

"""
    gertler_karadi_2015(; ndraws=200, chol_src=nothing, iv_src=nothing)

Gertler & Karadi (2015, AEJ Macro): a four-variable monthly VAR in the
1-year rate, log CPI, log industrial production and the excess bond premium,
12 lags, 48-step horizon, with wild-bootstrap bands.

Computes both the recursive benchmark and the external-instrument
identification, the latter using the FF4 surprise as the instrument.

!!! note
    `ident = "iv"` instruments the **first** column of `ENDO`, so `gs1` must
    be ordered first — as it is here, and as `GO_GK2015.m` orders it.
"""
function gertler_karadi_2015(; ndraws::Int = 200, chol_src = nothing,
                             iv_src = nothing, rng = Random.default_rng())
    data, vnames, _, _ = load_replication_data(joinpath("GK2015", "GK2015_Data.xlsx"))

    ENDO = hcat((col(data, vnames, n) for n in ("gs1", "logcpi", "logip", "ebp"))...)
    IV = reshape(col(data, vnames, "ff4_tc"), :, 1)

    VAR, VARopt = VARmodel(ENDO, 12, 1)
    VAR = @set VAR.IV = IV
    VARopt = @set VARopt.nsteps = 48
    VARopt = @set VARopt.ndraws = ndraws
    VARopt = @set VARopt.method = "wild"

    ochol = @set VARopt.ident = "short"
    IR_chol, VARc = VARir(VAR, ochol)
    cINF, cSUP, cMED, cBAR = VARirband(VARc, ochol; rng = rng, src = chol_src)

    oiv = @set VARopt.ident = "iv"
    IR_iv, VARi = VARir(VAR, oiv)
    iINF, iSUP, iMED, iBAR = VARirband(VARi, oiv; rng = rng, src = iv_src)

    return (; IR_chol, B_chol = VARc.B, cINF, cSUP, cMED, cBAR,
            IR_iv, B_iv = VARi.B, Biv = VARi.Biv, iINF, iSUP, iMED, iBAR,
            VAR, VARopt)
end

"""
    uhlig_2005(; ndraws=100, src=nothing, post_draws=nothing)

Uhlig (2005, JME): a six-variable monthly VAR identified by sign
restrictions on the monetary policy shock — the deflator, commodity prices
and non-borrowed reserves fall, the federal funds rate rises, and output and
total reserves are left unrestricted. The restrictions are imposed for six
months (`sr_hor = 6`) and the bands are 68% credible intervals (`pctg = 68`),
as in `GO_Uhlig2005.m`.

12 lags, 60-step horizon. `GO_Uhlig2005.m` rescales five of the six series
by 100 before estimation; that is reproduced here.

Defaults to `ndraws = 100` rather than the script's 500 purely to keep test
runtime and fixture size down.
"""
function uhlig_2005(; ndraws::Int = 100, src = nothing, post_draws = nothing,
                    rng = Random.default_rng())
    data, vnames, _, _ = load_replication_data(joinpath("Uhlig2005", "Uhlig2005_Data.xlsx"))

    # GO_Uhlig2005.m scales these five by 100 (the fed funds rate is left alone).
    for n in ("y", "pi", "comm", "nbres", "res")
        j = findfirst(==(n), vnames)
        j === nothing && continue
        data[:, j] .*= 100
    end

    VAR, VARopt = VARmodel(data, 12, 1)
    VARopt = @set VARopt.nsteps = 60
    VARopt = @set VARopt.ndraws = ndraws
    VARopt = @set VARopt.sr_hor = 6
    VARopt = @set VARopt.pctg = 68.0

    #        y   pi  comm nbres res  ff
    SIGN = [0.0 0 0 0 0 0     # Real GDP        -- unrestricted
            -1.0 0 0 0 0 0    # Deflator        -- falls
            -1.0 0 0 0 0 0    # Commodity price -- falls
            0.0 0 0 0 0 0     # Total reserves  -- unrestricted
            -1.0 0 0 0 0 0    # Nonborrowed res -- falls
            1.0 0 0 0 0 0]    # Fed funds rate  -- rises

    out = SR(VAR, SIGN, VARopt; rng = rng, src = src, post_draws = post_draws)
    return (; SRout = out, SIGN, VAR, VARopt)
end
