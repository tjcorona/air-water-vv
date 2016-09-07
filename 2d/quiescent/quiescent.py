from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import MeshRefinement as mr
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 1., "Height of free surface above bottom"),
    # tank
    ("tank_x", 4., "Dimensions of the tank"),
    ("tank_y", 2., "Dimensions of the tank"),
    # waves
    # caisson
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("T", 50.,"Simulation time"),
    ("he", 0.02,"mesh size"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.33 ,"Target cfl"),
    ("nsave",  1,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel"),
    # numerical
    ("ns_forceStrongDirichlet", False, ""),
    ("backgroundDiffusionFactor", 0.01, ""),
    ("shockCapturingFactor", 0.25, ""),
    ("sc_uref", 1.0, ""),
    ("sc_beta", 1.0, ""),
    ("useVF", 1.0, ""),
    ("lag_shockCapturing", True, ""),
    ("epsFact_density", 3.0, ""),
    ("epsFact_redistance", 0.33, ""),
    ("epsFact_consrv_diffusion", 1.0, ""),
    ("redist_Newton", True, ""),
    ("ns_lag_subgridError", True, "")])




# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level


tank_dim = (opts.tank_x, opts.tank_y)


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank = st.Rectangle(domain, tank_dim)
tank.translate((tank_dim[0]/2., tank_dim[1]/2.))


tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['x+'].setNoSlip()


he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)

##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=False  # to use moving mesh
checkMass=False  # ???
applyCorrection=True  # ???
applyRedistancing=True  # for redist
freezeLevelSet=True  # 

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T or 100*period
nDTout = int(T*opts.nsave)
if opts.nsave != 0:
    dt_out =  (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = False
spaceOrder = 1
useHex     = False  # discretization. True: cube, False: simplex
useRBLES   = 0.0  # multiplied with subgridError
useMetrics = 1.0  # ???
useVF = opts.useVF  # use in smoothing
useOnlyVF = False  # use of ls and rd or not
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988

#  Discretization
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters

ns_forceStrongDirichlet = opts.ns_forceStrongDirichlet
backgroundDiffusionFactor = opts.backgroundDiffusionFactor

ns_shockCapturingFactor = ls_shockCapturingFactor = rd_shockCapturingFactor = vof_shockCapturingFactor = opts.shockCapturingFactor
ls_sc_uref = vof_sc_uref = opts.sc_uref
ls_sc_beta = vof_sc_beta = opts.sc_beta
ns_lag_shockCapturing = rd_lag_shockCapturing = vof_lag_shockCapturing = ls_lag_shockCapturing = opts.lag_shockCapturing
epsFact_density    = opts.epsFact_density
epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density

epsFact_redistance = opts.epsFact_redistance
epsFact_consrv_diffusion = opts.epsFact_consrv_diffusion
redist_Newton = opts.redist_Newton
ns_lag_subgridError = opts.ns_lag_subgridError


# tolerances

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*he**2)
rd_nl_atol_res = max(1.0e-12,0.01*he)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank.dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))


model_n = -1
if movingDomain:
    MD_model = model_n+1  # index moving mesh model
    model_n += 1
else:
    MD_model = None

# Navier-Stokes and VOF
V_model = model_n+1  # index twp model
VF_model = model_n+2  # index vof model
model_n += 2
if not useOnlyVF:
    LS_model = model_n+1  # index ld model
    RD_model = model_n+2  # index rd model
    LSC_model = model_n+3  # index ls_consrv model
    model_n += 3
else:
    LS_model = None
    RD_model = None
    LSC_model = None


if useRANS > 0:
    K_model = model_n+1  # index kappa model
    D_model = model_n+2  # index dissipation model
    model_n += 2
