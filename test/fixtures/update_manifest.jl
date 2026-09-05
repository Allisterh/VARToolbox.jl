# =======================================================================
# update_manifest.jl
#
# Recompute test/fixtures/MANIFEST.toml.
#
# The manifest records a SHA-256 for every committed .mat fixture together
# with the generator that produced it. Both the generators and their output
# are committed, so `] test` runs on a machine with no MATLAB licence, while
# a fixture that has drifted from its generator is still detectable.
#
# Run after regenerating any fixture:
#
#   julia --project=test test/fixtures/update_manifest.jl
#
# Regenerating the fixtures themselves needs MATLAB:
#
#   matlab -batch "addpath('test/fixtures/gen'); gen_phase1('test/fixtures')"
# =======================================================================

using SHA
using Dates
using MAT
include(joinpath(@__DIR__, "content_hash.jl"))

const FIXDIR = @__DIR__

# fixture stem => generator that produces it
const GENERATORS = Dict(
    "phase1" => "gen_phase1.m",
    "phase2" => "gen_phase2.m",
    "phase3" => "gen_phase3.m",
    "phase3b" => "gen_phase3b.m",
    "phase4" => "gen_phase4.m",
    "phase5" => "gen_phase5.m",
    "phase6" => "gen_phase6.m",
    "phase7" => "gen_phase7.m",
    "phase9" => "gen_phase9.m",
)

function sha256_file(path)
    open(path) do io
        bytes2hex(sha256(io))
    end
end

function main()
    mats = sort(filter(f -> endswith(f, ".mat"), readdir(FIXDIR)))
    isempty(mats) && error("no .mat fixtures found in $FIXDIR")

    io = IOBuffer()
    println(io, "# Golden-master fixture manifest -- GENERATED, do not hand-edit.")
    println(io, "# Regenerate with: julia --project=test test/fixtures/update_manifest.jl")
    println(io, "#")
    println(io, "# Fixtures are committed alongside their generators so that the test")
    println(io, "# suite runs without a MATLAB licence.")
    println(io, "#")
    println(io, "# Two hashes per fixture, and they do different jobs:")
    println(io, "#")
    println(io, "#   content_sha256  a canonical hash over the DATA. Stable across")
    println(io, "#                   regeneration, so a generator that has drifted from")
    println(io, "#                   its committed output changes it. This is what the")
    println(io, "#                   test suite asserts.")
    println(io, "#   sha256          a hash of the file BYTES. Recorded for integrity,")
    println(io, "#                   but MATLAB stamps a creation time into every .mat")
    println(io, "#                   header, so this changes on EVERY regeneration even")
    println(io, "#                   when the data is bit-identical. Checked as a")
    println(io, "#                   warning only.")
    println(io)
    println(io, "[meta]")
    println(io, "updated = \"", Dates.format(Dates.now(), "yyyy-mm-dd"), "\"")
    println(io, "matlab = \"R2026a Update 3 (26.1.0.3276743)\"")
    println(io, "source_commit = \"7695feb9c72812e6f6f08c1288ccc8b927af0aa7\"")
    # phase9.mat alone is generated from VAR Toolbox 4.0 (fork dir v4dot0/):
    println(io, "source_commit_v4dot0 = \"36e706a0b6af6a6a901ce7084e5502ab71b75f42\"")
    println(io)

    for m in mats
        stem = splitext(m)[1]
        gen = get(GENERATORS, stem, "unknown")
        println(io, "[fixtures.", stem, "]")
        println(io, "file = \"", m, "\"")
        println(io, "generator = \"gen/", gen, "\"")
        println(io, "sha256 = \"", sha256_file(joinpath(FIXDIR, m)), "\"")
        println(io, "content_sha256 = \"", content_hash(matread(joinpath(FIXDIR, m))), "\"")
        println(io, "bytes = ", filesize(joinpath(FIXDIR, m)))
        println(io)
    end

    write(joinpath(FIXDIR, "MANIFEST.toml"), String(take!(io)))
    println("wrote MANIFEST.toml for ", length(mats), " fixture(s): ", join(mats, ", "))
end

main()
