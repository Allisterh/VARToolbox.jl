#!/usr/bin/env python3
"""Regenerate the fixture shims from the upstream sources.

Each shim is byte-identical to its upstream file except that its random
number call sites are routed through a recording wrapper, so the exact draw
sequence MATLAB consumed can be replayed in Julia (M2):

    SignRestrictions.m   randn( -> vartb_recrandn(     (2 sites)
    VARirband.m          rand(  -> vartb_recrand(      (3 sites)
    VARvdband.m          rand(  -> vartb_recrand(      (2 sites)

Byte mode throughout: the upstream files use CRLF line endings, and reading
them as text would silently rewrite them, making each shim look wholly
changed in a diff and destroying the audit trail that these are minimal edits.

Run:  python3 julia/VARToolbox/test/fixtures/gen/shim/make_shim.py
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[5]                      # shim/gen/fixtures/test/VARToolbox/julia -> repo root

# (upstream file, token to replace, replacement, expected number of sites)
JOBS = [
    ("SignRestrictions.m", b"randn(", b"vartb_recrandn(", 2),
    ("VARirband.m",        b"rand(",  b"vartb_recrand(",  3),
    ("VARvdband.m",        b"rand(",  b"vartb_recrand(",  2),
]

TEMPLATE = (
    b"%% ==== FIXTURE SHIM -- do not use outside fixture generation ====\r\n"
    b"%% Byte-identical to v3dot1/VAR/%s except that its %d\r\n"
    b"%% random-number call sites are routed through a recording wrapper, so\r\n"
    b"%% the exact draw sequence can be replayed in Julia (M2). Nothing else\r\n"
    b"%% is changed: algorithm, branch structure and acceptance logic are\r\n"
    b"%% untouched.\r\n"
    b"%% Regenerate with test/fixtures/gen/shim/make_shim.py.\r\n"
)

for name, token, repl, expected in JOBS:
    src_path = ROOT / "v3dot1" / "VAR" / name
    src = src_path.read_bytes()
    # randn( also contains rand( -- count only standalone occurrences.
    n = src.count(token) - (src.count(b"randn(") if token == b"rand(" else 0)
    if n != expected:
        raise SystemExit(f"{name}: expected {expected} {token.decode()} sites, found {n}")
    if token == b"rand(":
        # Protect randn( before rewriting rand(, then restore it.
        out = src.replace(b"randn(", b"\x00RANDN\x00")
        out = out.replace(token, repl).replace(b"\x00RANDN\x00", b"randn(")
    else:
        out = src.replace(token, repl)
    header = TEMPLATE % (name.encode(), expected)
    (HERE / name).write_bytes(header + out)
    print(f"wrote {HERE / name} ({n} sites routed)")

# VARdrawpost is different in kind: its randomness comes from wishrnd and
# mvnrnd, whose streams cannot be reproduced across languages at all. So
# rather than routing a call site, we append three lines that RECORD the
# draw it produced. The computation itself is untouched.
POST_HEADER = (
    b"%% ==== FIXTURE SHIM -- do not use outside fixture generation ====\r\n"
    b"%% Byte-identical to v3dot1/VAR/VARdrawpost.m with three lines appended\r\n"
    b"%% that record the draw it produced, so the identical sequence can be\r\n"
    b"%% injected into the Julia port (M2). The draw itself uses wishrnd and\r\n"
    b"%% mvnrnd, whose streams cannot be reproduced across languages -- so the\r\n"
    b"%% draws are recorded rather than replayed from a common seed.\r\n"
    b"%% Nothing in the computation is changed.\r\n"
    b"%% Regenerate with test/fixtures/gen/shim/make_shim.py.\r\n"
)
POST_TAIL = (
    b"\r\n"
    b"%% ---- fixture recording (appended by make_shim.py) ----\r\n"
    b"global VARTB_POST\r\n"
    b"if isempty(VARTB_POST); VARTB_POST = {}; end\r\n"
    b"VARTB_POST{end+1} = {sigma_draw, Ft_draw};\r\n"
)
_src = (ROOT / "v3dot1" / "VAR" / "VARdrawpost.m").read_bytes()
(HERE / "VARdrawpost.m").write_bytes(POST_HEADER + _src + POST_TAIL)
print(f"wrote {HERE / 'VARdrawpost.m'} (3 recording lines appended)")
