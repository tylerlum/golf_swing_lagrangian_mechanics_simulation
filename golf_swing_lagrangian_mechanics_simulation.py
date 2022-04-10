# -*- coding: utf-8 -*-
# # Golf Swing Lagrangian Mechanics Simulation
#
# The purpose of this notebook is to run numerical simulations of a golf swing modeled with Lagrangian mechanics.

import numpy as np
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
from enum import Enum


# ## Important Variables, Abbreviations, and Units
#
# We will be using the following abbreviations throughout this notebook:
#
# * L\_ = length
#
# * M\_ = mass
#
# * D\_ = derivative wrt time
#
# * \_A = arm
#
# * \_S = golf club shaft
#
# * \_C = golf clubhead
#
# Unless otherwise noted, we will be using SI units for all values in the simulation. For example:
#
# * angle \[radians\]
#
# * time \[seconds\]
#
# * length \[meters\]
#
# * mass \[kg\]
#
# * energy \[joules\]

# ## Select simulation parameters

class Forward_Sim_Method(Enum):
    FORWARD_EULER = 0
    IMPROVED_EULER = 1


DT = 1e-2
FORWARD_METHOD = Forward_Sim_Method.FORWARD_EULER


# ## Select golf swing variable parameters

class System_Type(Enum):
    PASSIVE_ARMS = 0
    CONTROLLED_ARMS = 1


# +
T_FIXED = 1  # Fixed wrist angle  for t in [0, t_fixed], then free wrist angle for t in [t_fixed, infty]
θ_0 = np.radians(90)  # Initial arm angle for backswing
SYSTEM_TYPE = System_Type.PASSIVE_ARMS
T_HORIZON = 10

M_A = 2 * 5.4  # Each arm about 5.4kg https://whatthingsweigh.com/how-much-does-an-arm-weigh/
L_A = 0.50 # Arm-length 0.5m https://www.craftyarncouncil.com/standards/man-size
M_S = 0.065 + 0.050  # Shaft weight about 65g + 50g https://www.hirekogolf.com/head-weights-shaft-weights-and-balance-points-oh-my
L_S = 45 * 2.54 / 100  # Shaft length about 45 inches https://www.hirekogolf.com/head-weights-shaft-weights-and-balance-points-oh-my
M_C = 0.2  # Clubhead about 200g https://www.hirekogolf.com/head-weights-shaft-weights-and-balance-points-oh-my
# -

# ## Setup fixed golf swing variable parameters

D_θ_0 = 0  # No initial arm angle speed at start of backswing
φ_0 = np.radians(90)  # Initial wrist angle is 90 degrees
D_φ_0 = 0  # No initial wrist angle speed at start of backswing

# ## Prepare variables for simulation

# +
n_steps = math.ceil(T_HORIZON / DT)

θ = np.zeros(n_steps)
D_θ = np.zeros(n_steps)
θ[0] = θ_0
D_θ[0] = D_θ_0
# -

# ## Run simulation

# +
for n in tqdm(range(1, n_steps)):
    # Compute all variables for this given time step
    t = n * DT
    
    
# -

# ## Visualize golf swing


