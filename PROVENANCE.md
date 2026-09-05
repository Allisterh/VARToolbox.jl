# Provenance

`VARToolbox.jl` is a Julia port of **Ambrogio Cesa-Bianchi's VAR Toolbox**
for MATLAB. It is a derivative work, not a clean-room reimplementation:
the algorithms, the option semantics, the array conventions, and in many
places the statement ordering are taken directly from the MATLAB source.

## Upstream

| | |
|---|---|
| Original author | Ambrogio Cesa-Bianchi (`ambrogiocesabianchi@gmail.com`) |
| Upstream repository | <https://github.com/ambropo/VAR-Toolbox> |
| Fork this port was made from | <https://github.com/hirothreading/VAR-Toolbox> |
| Source tree | `v3dot1/` — "VAR Toolbox 3.1" |
| Derived from commit | `7695feb9c72812e6f6f08c1288ccc8b927af0aa7` |
| Source size at that commit | 115 `.m` files, 11,417 lines |

### Upstream v4.0 re-sync (September 2026)

Upstream released **VAR Toolbox 4.0** on 2026-06-27 (`0ed54eaa1c1a…`) and
a follow-up bug-fix commit on 2026-08-18 (`36e706a0b6af6a6a901ce7084e5502ab71b75f42`),
rewriting its git history and flattening the tree in the process. The fork
carries that tree as `v4dot0/` (a git subtree of `upstream/main`). v4.0 is
an API redesign and is **not** what this package ports; but its
`ChangeLog.md` lists numerous bug fixes to computations that v3.1 and this
port share, and those fixes have been adopted:

| | |
|---|---|
| Upstream commit whose fixes are adopted | `36e706a0b6af6a6a901ce7084e5502ab71b75f42` |
| Fork directory | `v4dot0/` |
| Record of what was and was not adopted | `PORTING_LOG.md`, "Phase 9" |
| Fixture generated from v4.0 | `test/fixtures/phase9.mat` (`gen/gen_phase9.m`) |

Every other fixture is still generated from `v3dot1/` at the commit above.

Any future re-sync should `git fetch upstream` in the fork, merge with
`-X subtree=v4dot0`, read the new `v4dot0/ChangeLog.md` entries, and work
the delta through `PORTING_LOG.md` as Phase 9 did.

## Licence

The original toolbox is licensed **GNU GPL-3.0**. This port therefore is,
and must remain, GPL-3.0. See `LICENSE` at the root of this repository.

Every ported source file carries a header naming the MATLAB function it
derives from, so that the correspondence survives refactoring.

## Reference environment

Golden-master fixtures under `test/fixtures/` were generated with:

| | |
|---|---|
| MATLAB | R2026a Update 3 (26.1.0.3276743) |
| Statistics and Machine Learning Toolbox | 26.1 (R2026a) |
| Platform | macOS 26.6, Apple silicon |
| Julia (validation side) | 1.12.5 |

`test/fixtures/MANIFEST.toml` records a SHA-256 for every fixture along
with the generator that produced it, so drift between a generator and its
committed output is detectable without a MATLAB licence.

Running the test suite and examples requires neither MATLAB nor the
MATLAB fork — the replication `.xlsx` sheets are vendored under `data/`
(see `data/README.md`). Two maintenance tasks do require a checkout of the
fork above alongside MATLAB itself: regenerating fixtures
(`test/fixtures/gen/gen_phase1..7.m` read its `v3dot1/` directly;
`gen_phase9.m` reads `v4dot0/`) and the vendored-data drift check in
`runtests.jl`, which activates when the `VARTOOLBOX_V3DOT1` environment
variable points at a `v3dot1/` checkout.
