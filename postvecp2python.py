# =================================================================== 
#                        POSTVECP to PYTHON
# code written by Clement Baruteau (CB) and Lorenzo Valdettaro (LV).
# The code is adapted from the GUI code written by Vincent Prat and
# inspired from an IDL code originally written by Michel Rieutord and
# Lorenzo Valdettaro
# ===================================================================

# =========================================
#            TO DO LIST
# =========================================
# - calculation of Lyapunov coefficients?
# - plot shell-integrated dissipation for a range of directories
# =========================================

# =====================
# IMPORT GLOBAL VARIABLES
# =====================
import par

# =====================
# DISPLAY MERIDIONAL CUT
# =====================
if par.plot_zcut == 'Yes':
    from plot_zcut import *
    plotzcut()

# =====================
# DISPLAY MERIDIONAL CUT in 3D via PARAVIEW
# =====================
if par.plot_zcut_3D == 'Yes':
    from plot_zcut_3D import *
    plotzcut3D()

# =====================
# DISPLAY MODE's SPECTRAL CONTENT
# =====================
if par.plot_spectrum == 'Yes':
    from plot_spectrum import *
    plotspectrum()

# =====================
# DISPLAY EIGENFREQUENCIES (QZ RUN ONLY)
# =====================
if par.plot_qz_eigenfq == 'Yes':
    from plot_qz_eigenfq import *
    plotqzeigenfq()

# =====================
# DISPLAY TOTAL DISSIPATION vs. FORCING FREQUENCY
# =====================
if par.plot_total_dissipation == 'Yes':
    from plot_total_dissipation import *
    plottotaldissipation()
