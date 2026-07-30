# =======================================================================
# VARToolbox.jl
#
# A Julia port of Ambrogio Cesa-Bianchi's VAR Toolbox for MATLAB (v3.1).
#
# This is a derivative work of GPL-3.0 licensed code and is itself
# GPL-3.0. See PROVENANCE.md for the upstream commit this was ported from
# and the reference environment used to generate golden-master fixtures.
#
# Design decisions that shape the whole package (see ../DESIGN.md):
#
#   D1  Functions keep their MATLAB names verbatim (`VARmodel`, `VARir`,
#       `SR`, ...). Types take distinct suffixed names (`VARModel`,
#       `VAROptions`, ...) so that the MATLAB idiom `VAR = VARmodel(...)`
#       does not shadow a type binding.
#
#   D2  Where the MATLAB source is odd, the oddity is reproduced and pinned
#       with a test rather than silently improved. Every such case carries a
#       docstring admonition and a PORTING_LOG entry.
#
#   D3  `VAROptions` and `LPOptions` are IMMUTABLE. MATLAB passes structs by
#       value, and the toolbox relies on that: `SR.m` sets `VARopt.ident =
#       'sign'` and `SignRestrictions.m` sets `VARopt.nsteps =
#       VARopt.sr_hor` on their local copies. A mutable port would leak
#       those back to the caller and corrupt every downstream band and plot.
#       Update options with Setfield's `@set`.
# =======================================================================

module VARToolbox

using LinearAlgebra
using Random
using Statistics
using StatsBase
using Distributions
using Setfield
using Printf

# ---- Phase 1: foundations ---------------------------------------------
include("auxiliary.jl")
include("utils.jl")

# ---- Phase 2: estimation ----------------------------------------------
include("ols.jl")
include("var.jl")
include("ardl.jl")

# ---- Phase 3: identification ------------------------------------------
include("identification.jl")
include("signrestrictions.jl")

# ---- Phase 4: error bands --------------------------------------------
include("bands.jl")

# Auxiliary
export tdis_prb, tdis_inv, fdis_prb, beta_pdf, beta_inv
export nanmean, nansum, nanstd, nanvar, nansem, nanmedian, nanmax, nanmin

# Utils
export CommonSample, Num2NaN, NaN2Num, getqr, OrthNorm
export trimr, L, seqa, ptrend
export VARmakelags, VARmakexy

# Estimation
export OLSmodel, OLSOutput
export VARmodel, VARoption, VARModel, VAROptions, EquationResult
export ARDLmodel, ARDLModel

# Identification
export VARir, VARvd, VARhd, HDOutput
export SR, SROutput, SignRestrictions, VARdrawpost, DrawSource

# Error bands
export VARirband, VARvdband

# Re-exported so `@set VARopt.nsteps = 40` works without an extra `using`.
export @set

end # module
