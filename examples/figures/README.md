# Replication figures

Every figure the MATLAB toolbox produces for the four published
replications, regenerated from the Julia port by `../figures.jl`:

```
julia --project=julia/VARToolbox/examples -e '
    using Plots
    include("julia/VARToolbox/examples/figures.jl")
    generate_all_figures()'
```

The bootstrap and sign-restriction draws are replayed from the committed
golden-master fixture `test/fixtures/phase6.mat` (recorded from MATLAB
R2026a), so every curve below is numerically identical — not just similar —
to a MATLAB run with those draws. The rendering itself was compared
side-by-side against figures produced by MATLAB R2026a's own plot layer;
the accepted rendering differences are listed in `PORTING_LOG.md`
("Phase 8").

| file | upstream counterpart | produced by |
|---|---|---|
| `SW2001_IR_1..3.pdf` | `v3dot1/Replic/SW2001/IR_1..3.pdf` | `VARirplot` |
| `SW2001_VD.pdf` | `v3dot1/Replic/SW2001/VD.pdf` | `VARvdplot` |
| `BQ1989_IR_1..2.pdf` | `v3dot1/Replic/BQ1989/IR_1..2.pdf` | `VARirplot` |
| `BQ_Replication.pdf` | `v3dot1/Replic/BQ1989/BQ_Replication.pdf` | hand-rolled, `GO_BQ1989.m:65-88` |
| `Uhligh_Replication.pdf` | `v3dot1/Replic/Uhlig2005/Uhligh_Replication.pdf` | hand-rolled, `GO_Uhlig2005.m:91-102` (the filename typo is upstream's) |
| `GK_Replication.pdf` | `v3dot1/Replic/GK2015/GK_Replication.pdf` | hand-rolled, `GO_GK2015.m:79-118` |
| `SW_HD_1..3.pdf` | `v3dot1/Primer/graphics/SW_HD_1..3.pdf` | `VARhdplot` (Primer §9.6 — the only upstream exercise of the HD plotter and date axes) |

Two content-level caveats when comparing against the *committed upstream*
PDFs (as opposed to a fresh MATLAB run):

- **Uhlig runs 100 accepted draws, not `GO_Uhlig2005.m`'s 500** — the
  committed fixture caps the recorded-draw count for repository size
  (`test/fixtures/gen/gen_phase6.m`). Bands are slightly less smooth; the
  identification content is unchanged.
- **The upstream PDFs embed different random draws** (they were rendered by
  the original author years ago), so their error bands wiggle differently.
  The exact-match comparison target is a fresh MATLAB run replaying the
  fixture draws, which is what the Phase 8 validation used.
