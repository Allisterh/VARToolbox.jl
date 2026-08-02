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

Any future re-sync against upstream should diff `v3dot1/` at the commit
above against the then-current upstream `v3dot1/`, and work the resulting
delta through `PORTING_LOG.md`.

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
monorepo — the replication `.xlsx` sheets are vendored under `data/`
(see `data/README.md`). Two maintenance tasks do require a checkout of the
monorepo fork above alongside MATLAB: regenerating fixtures
(`test/fixtures/gen/*.m` read `v3dot1/` directly) and the vendored-data
drift check in `runtests.jl`, which activates automatically whenever
`v3dot1/` is present three levels up.
