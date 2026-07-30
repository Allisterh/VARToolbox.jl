#!/usr/bin/env python3
"""Regenerate the SignRestrictions fixture shim from the upstream source.

The shim is byte-identical to v3dot1/VAR/SignRestrictions.m except that its
two randn() call sites are routed through vartb_recrandn(), which records
each draw so the identical sequence can be replayed in Julia (M2).

Byte mode throughout: the upstream file uses CRLF line endings, and reading
it as text would silently rewrite them, making the shim look wholly changed
in a diff and destroying the audit trail that this is a two-line edit.

Run:  python3 julia/VARToolbox/test/fixtures/gen/shim/make_shim.py
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[5]                      # shim/gen/fixtures/test/VARToolbox/julia -> repo root
SRC = ROOT / "v3dot1" / "VAR" / "SignRestrictions.m"
DST = HERE / "SignRestrictions.m"

HEADER = (
    b"% ==== FIXTURE SHIM -- do not use outside fixture generation ====\r\n"
    b"% Byte-identical to v3dot1/VAR/SignRestrictions.m except that its two\r\n"
    b"% randn() call sites are routed through vartb_recrandn(), which records\r\n"
    b"% each draw so the identical sequence can be replayed in Julia (M2).\r\n"
    b"% Nothing else is changed: the algorithm, the branch structure and the\r\n"
    b"% acceptance logic are untouched.\r\n"
    b"% Regenerate with test/fixtures/gen/shim/make_shim.py if upstream changes.\r\n"
)

src = SRC.read_bytes()
n = src.count(b"randn(")
if n != 2:
    raise SystemExit(f"expected 2 randn call sites in {SRC}, found {n}")
DST.write_bytes(HEADER + src.replace(b"randn(", b"vartb_recrandn("))
print(f"wrote {DST} ({n} randn sites routed)")
