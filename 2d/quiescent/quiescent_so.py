"""
Split operator module for two-phase flow
"""

import os
from proteus.default_so import *
from proteus import Context

# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise NameError, 'Split operator module must end with "_so.py"'

try:
    case = __import__(name)
    Context.setFromModule(case)
    ct = Context.get()
except ImportError:
    raise ImportError, str(name) + '.py not found'

for bc in ct.domain.bc:
    bc.getContext()
    
# List of p/n files
pnList = []
model_n = -1

# moving mesh
if ct.movingDomain:
    pnList += [("moveMesh_p", "moveMesh_n")]
    MD_model = model_n+1  # index moving mesh model
    model_n += 1
else:
    MD_model = None

# Navier-Stokes and VOF
pnList += [("twp_navier_stokes_p", "twp_navier_stokes_n"),
           ("vof_p", "vof_n")]
V_model = model_n+1  # index twp model
VF_model = model_n+2  # index vof model
model_n += 2

# Level set
if not ct.useOnlyVF:
    pnList += [("ls_p", "ls_n"),
               ("redist_p", "redist_n"),
               ("ls_consrv_p", "ls_consrv_n")]
    LS_model = model_n+1  # index ld model
    RD_model = model_n+2  # index rd model
    LSC_model = model_n+3  # index ls_consrv model
    model_n += 3
else:
    LS_model = None
    RD_model = None
    LSC_model = None

# Turbulence
if ct.useRANS > 0:
    pnList += [("kappa_p", "kappa_n"),
               ("dissipation_p", "dissipation_n")]
    K_model = model_n+1  # index kappa model
    D_model = model_n+2  # index dissipation model
    model_n += 2

#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

if ct.opts.nsave == 0:
    archiveFlag = ArchiveFlags.EVERY_USER_STEP
    tnList = [0., ct.dt_init, ct.T]
else:
    tnList=[0.0,ct.dt_init]+[ct.dt_init+ i*ct.dt_out for i in range(1,ct.nDTout+1)]
