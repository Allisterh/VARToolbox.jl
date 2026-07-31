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
    println(io, "# suite runs without a MATLAB licence. These hashes make drift between")
    println(io, "# a generator and its committed output detectable.")
    println(io)
    println(io, "[meta]")
    println(io, "updated = \"", Dates.format(Dates.now(), "yyyy-mm-dd"), "\"")
    println(io, "matlab = \"R2026a Update 3 (26.1.0.3276743)\"")
    println(io, "source_commit = \"7695feb9c72812e6f6f08c1288ccc8b927af0aa7\"")
    println(io)

    for m in mats
        stem = splitext(m)[1]
        gen = get(GENERATORS, stem, "unknown")
        println(io, "[fixtures.", stem, "]")
        println(io, "file = \"", m, "\"")
        println(io, "generator = \"gen/", gen, "\"")
        println(io, "sha256 = \"", sha256_file(joinpath(FIXDIR, m)), "\"")
        println(io, "bytes = ", filesize(joinpath(FIXDIR, m)))
        println(io)
    end

    write(joinpath(FIXDIR, "MANIFEST.toml"), String(take!(io)))
    println("wrote MANIFEST.toml for ", length(mats), " fixture(s): ", join(mats, ", "))
end

main()
