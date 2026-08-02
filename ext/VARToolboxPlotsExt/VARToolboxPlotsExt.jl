# =======================================================================
# VARToolboxPlotsExt — the plotting layer of VARToolbox.jl
#
# Package extension, activated by `using Plots`. Ports the figure-producing
# surface of VAR Toolbox 3.1 (Ambrogio Cesa-Bianchi): the six plot functions
# in VAR/ (VARirplot, VARvdplot, VARhdplot, SRirplot, SRvdplot, SRhdplot)
# and the Figure/ primitives they stand on (PlotSwathe, AreaPlot, cmap,
# DatesPlot, FigSize, SaveFigure).
#
# Fidelity target: v3dot1 logic with v3dot0's `quality==2` save branches
# restored — v3dot1 dropped them while defaulting `quality = 2`, so the
# stock v3dot1 plot layer saves nothing. The toolbox's own committed
# reference figures embody the v3dot0 behaviour. Every deviation is listed
# in PORTING_LOG.md ("Phase 8").
#
# GPL-3.0, derivative of https://github.com/ambropo/VAR-Toolbox.
# =======================================================================

module VARToolboxPlotsExt

using VARToolbox
using Plots
using Plots: RGB
using VARToolbox: VAROptions, HDOutput
import VARToolbox: VARirplot, VARvdplot, VARhdplot, SRirplot, SRvdplot, SRhdplot
import VARToolbox: PlotSwathe, PlotSwathe!, AreaPlot!, cmap, DatesPlot!,
                   FigSize, SaveFigure

include("primitives.jl")
include("varplots.jl")

end # module
