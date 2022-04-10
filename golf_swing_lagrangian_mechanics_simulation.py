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
    FORWARD_EULER = "FORWARD_EULER"
    IMPROVED_EULER = "IMPROVED_EULER"


DT = 1e-2
FORWARD_METHOD = Forward_Sim_Method.FORWARD_EULER


# ## Select golf swing variable parameters

class Arm_Type(Enum):
    PASSIVE_ARMS = "PASSIVE_ARMS"
    CONTROLLED_ARMS = "CONTROLLED_ARMS"


class Wrist_Type(Enum):
    FIXED_WRIST = "FIXED_WRIST"
    LOOSE_WRIST = "LOOSE_WRIST"


# +
T_FIXED = 1  # Fixed wrist angle  for t in [0, t_fixed], then free wrist angle for t in [t_fixed, infty]
θ_0 = np.radians(90)  # Initial arm angle for backswing
ARM_TYPE = Arm_Type.PASSIVE_ARMS
T_HORIZON = 10
CONTROLLED_ARM_ANGULAR_ACCEL = 1

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
g = 9.8  # Acceleration due to gravity

# ## Prepare variables for simulation

# +
n_steps = math.ceil(T_HORIZON / DT)
t = np.arange(n_steps) * DT

# Variables to be populated
θ = np.zeros(n_steps)
D_θ = np.zeros(n_steps)
φ = np.zeros(n_steps)
D_φ = np.zeros(n_steps)

# Initial conditions
θ[0] = θ_0
D_θ[0] = D_θ_0
φ[0] = φ_0
D_φ[0] = D_φ_0


# -

# ## Arm and Wrist Specific Forward Simulation

# +
def passive_arm_fixed_wrist_dd_θ_f(θ, D_θ, φ, D_φ, t):
    return ((-(M_A/2 + M_S + M_C)*L_A + (M_S/2 + M_C)*L_S) * g * np.sin(θ) /
            ((M_A/3 + M_S + M_C)*L_A**2 + (M_S/3 + M_C)*L_S**2))
    
def passive_arm_fixed_wrist_dd_φ_f(θ, D_θ, φ, D_φ, t):
    return 0

def _passive_arm_loose_wrist_helper(θ, D_θ, φ, D_φ, t):
    # E*dd_θ + F*dd_φ = Z
    E = (M_A/3 + M_S + M_C)*L_A**2 + (M_S/3 + M_C)*L_S**2 - (M_S + 2*M_C)*L_A*L_S*np.cos(φ)
    F = -(M_S/3 + M_C)*L_S**2 + (M_S/2 + M_C)*L_A*L_S*np.cos(φ)
    G = (M_S + 2*M_C)*L_A*L_S*np.sin(φ)
    H = -(M_S/2 + M_C)*L_A*L_S_np.sin(φ)
    Z = -(M_A/2 + M_S + M_C)*g_L_A*np.sin(θ) - (M_S/2 + M_C)*g*L_S*np.sin(φ - θ) - G*D_θ*D_φ - H*D_φ**2
    
    # J*dd_θ + K*dd_φ = P
    J = (M_S/2 + M_C)*L_A*L_S*np.cos(φ) - (M_S/3 + M_C)*L_S**2
    K = (M_S/3 + M_C)*L_S**2
    P = -(M_S/2 + M_C)*L_A*L_S*D_θ*(D_φ - D_θ)*np.sin(φ) + (M_S/2 + M_C)*g*L_S*np.sin(φ - θ)
    
    # Solve system of equations
    # A = (E F
    #      J K)
    # x = (dd_θ
    #      dd_φ)
    # b = (Z
    #      P)
    A = np.array([[E, F],
                  [J, K]])
    b = np.array([Z,
                  P])
    x = np.squeeze(np.linalg.solve(A, b))
    assert(x.size == 2)
    
    dd_θ, dd_φ = x[0], x[1]
    return dd_θ, dd_φ

def passive_arm_loose_wrist_dd_θ_f(θ, D_θ, φ, D_φ, t):
    dd_θ, _ = _passive_arm_loose_wrist_helper(θ, D_θ, φ, D_φ, t)
    return dd_θ

def passive_arm_loose_wrist_dd_φ_f(θ, D_θ, φ, D_φ, t):
    _, dd_φ = _passive_arm_loose_wrist_helper(θ, D_θ, φ, D_φ, t)
    return dd_φ

def controlled_arm_fixed_wrist_dd_θ_f(θ, D_θ, φ, D_φ, t):
    return CONTROLLED_ARM_ANGULAR_ACCEL

def controlled_arm_fixed_wrist_dd_φ_f(θ, D_θ, φ, D_φ, t):
    return 0

def controlled_arm_loose_wrist_dd_θ_f(θ, D_θ, φ, D_φ, t):
    return CONTROLLED_ARM_ANGULAR_ACCEL

def controlled_arm_loose_wrist_dd_φ_f(θ, D_θ, φ, D_φ, t):
    dd_θ = controlled_arm_loose_wrist_dd_θ_f(θ, D_θ, φ, D_φ, t)
    W = -(M_S/2 + M_C)*L_A*L_S*D_θ*(D_φ - D_θ)*np.sin(φ) + (M_S/2 + M_C)*g*L_S*np.sin(φ - θ)
    numerator = W + (M_S/3 + M_C)*L_S**2*dd_θ - (M_S/2 + M_C)*L_A*L_S*dd_θ*np.cos(φ)
    denominator = (M_S/3 + M_C)*L_S**2
    return numerator / denominator


# -

# ## Run simulation

for n in tqdm(range(n_steps - 1)):
    # Currently at step n (time t), will be calculating values for n+1 (time t+dt)
    
    # Validate input
    wrist_type = Wrist_Type.FIXED_WRIST if t < T_FIXED else Wrist_Type.LOOSE_WRIST
    if wrist_type != Wrist_Type.FIXED_WRIST and wrist_type != Wrist_Type.LOOSE_WRIST:
        raise ValueError(f"Invalid wrist_type = {wrist_type}")
    if ARM_TYPE != Arm_Type.PASSIVE_ARMS and ARM_TYPE != Arm_Type.CONTROLLED_ARMS:
        raise ValueError(f"Invalid ARM_TYPE = {ARM_TYPE}")

    # Get arm type and wrist type specific functions
    if ARM_TYPE == Arm_Type.PASSIVE_ARMS and wrist_type == Wrist_Type.FIXED_WRIST:
        dd_θ_f = passive_arm_fixed_wrist_dd_θ_f
        dd_φ_f = passive_arm_fixed_wrist_dd_φ_f
    elif ARM_TYPE == Arm_Type.PASSIVE_ARMS and wrist_type == Wrist_Type.LOOSE_WRIST:
        dd_θ_f = passive_arm_loose_wrist_dd_θ_f
        dd_φ_f = passive_arm_loose_wrist_dd_φ_f
    elif ARM_TYPE == Arm_Type.CONTROLLED_ARMS and wrist_type == Wrist_Type.FIXED_WRIST:
        dd_θ_f = controlled_arm_fixed_wrist_dd_θ_f
        dd_φ_f = controlled_arm_fixed_wrist_dd_φ_f
    elif ARM_TYPE == Arm_Type.CONTROLLED_ARMS and wrist_type == Wrist_Type.LOOSE_WRIST:
        dd_θ_f = controlled_arm_loose_wrist_dd_θ_f
        dd_φ_f = controlled_arm_loose_wrist_dd_φ_f
    else:
        raise ValueError(f"Invalid pair ARM_TYPE = {ARM_TYPE}, wrist_type = {wrist_type}")

    # Calculate update
    if FORWARD_METHOD == Forward_Sim_Method.FORWARD_EULER:
        D_θ[n+1] = D_θ[n] + dt * dd_θ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n])
        θ[n+1] = θ[n] + dt * D_θ[n]
        D_φ[n+1] = D_φ[n] + dt * dd_φ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n])
        φ[n+1] = φ[n] + dt * D_φ[n]
    elif FORWARD_METHOD == Forward_Sim_Method.IMPROVED_EULER:
        raise NotImplementedError()
        # Compute intermediate values
        D_θ_n_plus_1 = D_θ[n] + dt * dd_θ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n])
        θ_n_plus_1 = θ[n] + dt * D_θ[n]
        D_φ_n_plus_1 = D_φ[n] + dt * dd_φ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n])
        φ_n_plus_1 = φ[n] + dt * D_φ[n]

        # Calculate update
        D_θ[n+1] = D_θ[n] + dt/2 * (dd_θ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n]) +
                                    dd_θ_f(θ_n_plus_1, D_θ_n_plus_1, φ_n_plus_1, D_φ_n_plus_1, t[n+1]))
        θ[n+1] = θ[n] + dt/2 * (D_θ[n] + D_θ[n+1])
        D_φ[n+1] = D_φ[n] + dt/2 * (dd_φ_f(θ[n], D_θ[n], φ[n], D_φ[n], t[n]) +
                                    dd_φ_f(θ_n_plus_1, D_θ_n_plus_1, φ_n_plus_1, D_φ_n_plus_1, t[n+1]))
        φ[n+1] = φ[n] + dt/2 * (D_φ[n] + D_φ[n+1])
    else:
        raise ValueError(f"Invalid FORWARD_METHOD = {FORWARD_METHOD}")

# ## Plot golf swing



# ## Visualize golf swing


