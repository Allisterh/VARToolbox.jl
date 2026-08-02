# Vendored replication data

Byte-for-byte copies of the `.xlsx` data files shipped inside Ambrogio
Cesa-Bianchi's VAR Toolbox (GPL-3.0, <https://github.com/ambropo/VAR-Toolbox>),
vendored so this package's examples and test suite run outside the
monorepo that also carries the MATLAB original:

| file | upstream source | used by |
|---|---|---|
| `SW2001/SW2001_Data.xlsx` | `v3dot1/Replic/SW2001/` | Stock & Watson (2001) replication |
| `BQ1989/BQ1989_Data.xlsx` | `v3dot1/Replic/BQ1989/` | Blanchard & Quah (1989) replication |
| `GK2015/GK2015_Data.xlsx` | `v3dot1/Replic/GK2015/` | Gertler & Karadi (2015) replication |
| `Uhlig2005/Uhlig2005_Data.xlsx` | `v3dot1/Replic/Uhlig2005/` | Uhlig (2005) replication |
| `Primer/SW2001_Data.xlsx` | `v3dot1/Primer/data/` | Primer §9 historical-decomposition figures (a *different* sheet from the Replic one — reordered variables, longer sample) |

Sheets use the toolbox's layout: long names on row 1, short names on row 2,
dates in column 1, and the sentinel `123456789` for missing values
(converted by `Num2NaN`).

Do not edit these. With the `VARTOOLBOX_V3DOT1` environment variable
pointing at a checkout of the MATLAB toolbox's `v3dot1/` directory,
`test/runtests.jl` asserts each copy is still byte-identical to its
upstream original, so drift fails the suite.
