using Test
using VARToolbox
using LinearAlgebra
using MAT
using SHA
include(joinpath(@__DIR__, "fixtures", "content_hash.jl"))
using Random
using Statistics
using Distributions

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
        # Fixtures are committed so the suite runs without MATLAB. Two
        # different checks, because they answer different questions:
        #
        #   content_sha256 -- has the DATA changed? Asserted. This is what
        #       catches a generator that drifted from its committed output.
        #   sha256 (bytes) -- has the FILE changed? Warned about only.
        #       MATLAB stamps a creation timestamp into every .mat header,
        #       so this differs after any regeneration even when every value
        #       is bit-identical (verified: regenerating phase1.mat gave a
        #       different file hash with all 115 fields comparing equal).
        manifest = joinpath(FIXDIR, "MANIFEST.toml")
        @test isfile(manifest)
        bytes_h = Dict{String,String}()
        content_h = Dict{String,String}()
        cur = ""
        for line in eachline(manifest)
            m = match(r"^\[fixtures\.(\w+)\]", line)
            m !== nothing && (cur = m.captures[1])
            isempty(cur) && continue
            b = match(r"^sha256 = \"([0-9a-f]{64})\"", line)
            b !== nothing && (bytes_h[cur] = b.captures[1])
            c = match(r"^content_sha256 = \"([0-9a-f]{64})\"", line)
            c !== nothing && (content_h[cur] = c.captures[1])
        end
        @test !isempty(content_h)
        @test keys(content_h) == keys(bytes_h)

        for (stem, want) in content_h
            path = joinpath(FIXDIR, "$stem.mat")
            @test isfile(path)
            got = content_hash(matread(path))
            if got != want
                @warn "fixture $stem.mat CONTENT differs from MANIFEST.toml -- " *
                      "the data changed. Regenerate with gen/gen_$stem.m and " *
                      "rerun update_manifest.jl, then review the diff."
            end
            @test got == want

            gotb = open(io -> bytes2hex(sha256(io)), path)
            gotb == bytes_h[stem] || @info(
                "fixture $stem.mat bytes differ from MANIFEST.toml but its " *
                "content matches -- expected after regeneration, since MATLAB " *
                "stamps a creation time into the .mat header.")
        end

        @testset "content hash is sensitive to the data" begin
            # A hash that never changes would pass every drift check
            # vacuously. Assert it actually moves when a value moves, when
            # a key is added or removed, and when NaN is involved -- but
            # NOT when only the file's byte layout differs.
            base = matread(joinpath(FIXDIR, "phase1.mat"))
            h0 = content_hash(base)

            perturbed = deepcopy(base)
            perturbed["A"][1, 1] += 1e-12
            @test content_hash(perturbed) != h0

            added = deepcopy(base)
            added["__extra__"] = [1.0]
            @test content_hash(added) != h0

            removed = deepcopy(base)
            delete!(removed, first(sort(collect(keys(removed)))))
            @test content_hash(removed) != h0

            # NaN must hash consistently -- fixtures are full of it.
            nanned = deepcopy(base)
            nanned["A"][2, 2] = NaN
            @test content_hash(nanned) != h0
            @test content_hash(nanned) == content_hash(deepcopy(nanned))

            # Key order must not matter.
            reordered = Dict{String,Any}(reverse(collect(base)))
            @test content_hash(reordered) == h0
        end
    end

    @testset "vendored data integrity" begin
        # data/ carries byte-for-byte copies of the toolbox's .xlsx sheets
        # so the package is self-contained outside the MATLAB monorepo
        # (https://github.com/hirothreading/VAR-Toolbox). To assert the
        # copies have not drifted from the originals, point
        # VARTOOLBOX_V3DOT1 at a checkout's v3dot1/ directory:
        #
        #   VARTOOLBOX_V3DOT1=~/GitHub/VAR-Toolbox/v3dot1 julia --project -e 'using Pkg; Pkg.test()'
        #
        # Without it the check skips (there is nothing to compare against).
        # The relative fallback covers the historical monorepo layout where
        # this package lived at <monorepo>/julia/VARToolbox.
        datadir = joinpath(@__DIR__, "..", "data")
        v3 = get(ENV, "VARTOOLBOX_V3DOT1",
                 normpath(joinpath(@__DIR__, "..", "..", "..", "v3dot1")))
        pairs = [
            ("SW2001", joinpath("Replic", "SW2001")),
            ("BQ1989", joinpath("Replic", "BQ1989")),
            ("GK2015", joinpath("Replic", "GK2015")),
            ("Uhlig2005", joinpath("Replic", "Uhlig2005")),
            ("Primer", joinpath("Primer", "data")),
        ]
        for (sub, upstream) in pairs
            vendored = readdir(joinpath(datadir, sub))
            files = filter(f -> endswith(f, ".xlsx"), vendored)
            @test !isempty(files)
            for f in files
                src = joinpath(v3, upstream, f)
                if isfile(src)
                    @test read(joinpath(datadir, sub, f)) == read(src)
                else
                    @info "MATLAB original not found (set VARTOOLBOX_V3DOT1); skipping drift check for data/$sub/$f"
                end
            end
        end
    end

    include("test_phase1.jl")
    include("test_phase2.jl")
    include("test_phase3.jl")
    include("test_phase3b.jl")
    include("test_phase4.jl")
    include("test_phase5.jl")
    include("test_phase6.jl")
    include("test_phase7.jl")
    include("test_phase8.jl")
    include("test_phase9.jl")
end
