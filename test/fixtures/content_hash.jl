# =======================================================================
# content_hash.jl
#
# A SHA-256 over the *data* in a .mat fixture, independent of the file's
# byte layout.
#
# Why this exists. MATLAB stamps every .mat file with a creation time in its
# 128-byte header:
#
#     MATLAB 5.0 MAT-file, Platform: MACA64, Created on: Fri Jul 31 00:20:32 2026
#
# so a byte-level hash of a regenerated fixture ALWAYS differs, even when
# every value in it is bit-identical. Verified: regenerating phase1.mat from
# an unchanged generator produced a different file hash while all 115 fields
# compared equal.
#
# That makes a byte hash useless for the thing the manifest is supposed to
# catch -- a generator that has drifted away from its committed output. It
# is still worth recording for file integrity, but the assertion the test
# suite makes has to be over content.
#
# The encoding below is canonical: keys sorted, types tagged, sizes
# included, and numeric data hashed from its raw bytes. NaN hashes
# consistently because MATLAB and Julia agree on the quiet-NaN bit pattern;
# `isequal` semantics (NaN == NaN) are what we want here and what the byte
# encoding gives.
# =======================================================================

using SHA

_hupdate!(ctx, s::AbstractString) = update!(ctx, Vector{UInt8}(codeunits(String(s))))
_hupdate!(ctx, n::Integer) = update!(ctx, reinterpret(UInt8, [Int64(n)]))

function _hash_value!(ctx, v)
    if v isa AbstractString
        _hupdate!(ctx, "str:")
        _hupdate!(ctx, v)
    elseif v isa Bool
        _hupdate!(ctx, "bool:")
        _hupdate!(ctx, Int(v))
    elseif v isa Number
        _hupdate!(ctx, "num:")
        update!(ctx, reinterpret(UInt8, [Float64(v)]))
    elseif v isa AbstractArray
        _hupdate!(ctx, "arr:")
        _hupdate!(ctx, ndims(v))
        for d in size(v)
            _hupdate!(ctx, d)
        end
        if eltype(v) <: Union{Float64,Float32,Int64,Int32,Bool}
            # Dense numeric block: hash the raw bytes of a Float64 copy so
            # that Int/Float storage differences across MATLAB versions do
            # not register as a content change. `vec` first -- a
            # comprehension over a matrix keeps its shape, and reinterpret
            # of a 2-D array is not the AbstractVector{UInt8} update! wants.
            flat = vec(Float64[Float64(x) for x in v])
            update!(ctx, collect(reinterpret(UInt8, flat)))
        else
            # Cell arrays and friends: recurse in column-major order.
            for x in v
                _hash_value!(ctx, x)
            end
        end
    elseif v isa AbstractDict
        _hupdate!(ctx, "dict:")
        for k in sort(collect(keys(v)))
            _hupdate!(ctx, string(k))
            _hupdate!(ctx, "=")
            _hash_value!(ctx, v[k])
        end
    elseif v === nothing || v === missing
        _hupdate!(ctx, "nil:")
    else
        _hupdate!(ctx, "other:")
        _hupdate!(ctx, string(typeof(v)))
        _hupdate!(ctx, string(v))
    end
    _hupdate!(ctx, ";")
    return ctx
end

"""
    content_hash(d::AbstractDict) -> String

Canonical SHA-256 over the data of a loaded `.mat` fixture. Stable across
regeneration with an unchanged generator, unlike a hash of the file bytes.
"""
function content_hash(d::AbstractDict)
    ctx = SHA2_256_CTX()
    _hash_value!(ctx, d)
    return bytes2hex(digest!(ctx))
end
