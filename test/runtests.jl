using Test
using VARToolbox
using LinearAlgebra
using MAT
using SHA

const FIXDIR = joinpath(@__DIR__, "fixtures")

"""
    fixture(name)

Load a golden-master fixture generated from the MATLAB original. See
`fixtures/gen/` for the generators and `PROVENANCE.md` for the reference
environment.
"""
fixture(name) = matread(joinpath(FIXDIR, "$name.mat"))

# Golden-master tolerance. Anything looser than this needs a written
# justification at the call site (see the porting brief's validation policy).
const RTOL = 1e-10

"MATLAB always returns matrices; collapse a fixture to a flat vector for comparison."
flat(x) = vec(collect(x))

"""
    approxeq(a, b; rtol=RTOL)

Elementwise comparison against a MATLAB fixture, with two deliberate
departures from plain `≈`.

**`NaN` matches `NaN`.** Plain `≈` returns false for `NaN ≈ NaN`, and NaN is
meaningful data throughout this toolbox: `CommonSample` counts it, `L` pads
with it, `HD` NaN-pads its first `nlag` rows, and the `nan*` family returns
it for all-NaN slices. Comparing those with `≈` reports spurious failures.

**An element passes if it is close either relatively or on the scale of the
array it belongs to**, i.e. `|x-y| <= rtol*max(|x|,|y|) + rtol*maxabs(b)`.

The absolute term matters. A VAR coefficient matrix routinely spans four
orders of magnitude — on the Stock-Watson data, coefficients run from
0.0042 to 1.97. A purely elementwise-relative test demands that the
smallest coefficient agree to `rtol` in its *own* units, which no two
linear-algebra implementations will do: MATLAB and Julia agree there to
1e-12 absolute, which is 1e-10 relative on that element but 1e-12 relative
to the matrix. That is a rounding artefact of the solve, not a porting
error.

A norm-based `isapprox(a, b; rtol)` would also pass, but it can hide a
single badly wrong element inside a large array. The scaled-atol form keeps
elementwise resolution while measuring "small" against the right yardstick:
a coefficient that is genuinely wrong by 1% fails, because its error is
large relative to the array maximum too.
"""
function approxeq(a, b; rtol = RTOL)
    va, vb = flat(a), flat(b)
    length(va) == length(vb) || return false
    scale = isempty(vb) ? 0.0 : maximum(x -> isfinite(x) ? abs(x) : 0.0, vb)
    atol = rtol * scale
    for (x, y) in zip(va, vb)
        if isnan(x) || isnan(y)
            (isnan(x) && isnan(y)) || return false
        elseif isinf(x) || isinf(y)
            x == y || return false
        else
            isapprox(x, y; rtol = rtol, atol = atol) || return false
        end
    end
    return true
end

macro golden(jl, mat)
    quote
        @test approxeq($(esc(jl)), $(esc(mat)))
    end
end

@testset "VARToolbox.jl" begin
    @testset "fixture integrity" begin
        # Fixtures are committed so the suite runs without MATLAB. Verify
        # they still hash to what the manifest recorded, so a fixture that
        # drifted from its generator fails loudly rather than silently
        # changing what "golden" means.
        manifest = joinpath(FIXDIR, "MANIFEST.toml")
        @test isfile(manifest)
        recorded = Dict{String,String}()
        cur = ""
        for line in eachline(manifest)
            m = match(r"^\[fixtures\.(\w+)\]", line)
            m !== nothing && (cur = m.captures[1])
            h = match(r"^sha256 = \"([0-9a-f]{64})\"", line)
            h !== nothing && !isempty(cur) && (recorded[cur] = h.captures[1])
        end
        @test !isempty(recorded)
        for (stem, want) in recorded
            path = joinpath(FIXDIR, "$stem.mat")
            @test isfile(path)
            got = open(path) do io
                bytes2hex(sha256(io))
            end
            if got != want
                @warn "fixture $stem.mat does not match MANIFEST.toml; " *
                      "regenerate with gen/gen_$stem.m then rerun update_manifest.jl"
            end
            @test got == want
        end
    end

    include("test_phase1.jl")
    include("test_phase2.jl")
end
