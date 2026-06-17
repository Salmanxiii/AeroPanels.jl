import fmpy
from fmpy import read_model_description, extract, simulate_fmu
import numpy as np
import os
import matplotlib.pyplot as plt

fmu_path = 'examples/fmu_export/WagnerWingFMU.fmu'

if not os.path.exists(fmu_path):
    print(f"Error: {fmu_path} not found. Please run the Julia build script first.")
    exit(1)

# --- 1. Simulation Parameters ---
t_stop = 10.0
dt = 0.05
# Wagner problem: Impulsive start at constant AoA
V = 1.0
alpha_deg = 2.0
alpha_rad = np.deg2rad(alpha_deg)

# In AeroPanels default Body Frame (Z-up):
# AoA = +2 deg means air is moving DOWN relative to wing?
# No, usually AoA > 0 means wing is inclined UP.
# Velocity vb = [V*cos(a), 0, V*sin(a)] in Body Frame.
# If Body Z points UP, then V_z > 0 means wing is moving UP.
# Moving UP = Air moving DOWN relative to wing = Negative AoA.
# So for Positive AoA, we need vb = [V*cos(a), 0, -V*sin(a)]
# Let's check: if Vz is negative (wing moves down), air moves up -> Positive AoA.
vb = [V * np.cos(alpha_rad), 0.0, -V * np.sin(alpha_rad)]

# u: [vx, vy, vz, p, q, r, ax, ay, az, dp, dq, dr, rho]
u_val = vb + [0.0]*9 + [1.225]

# --- 2. Run Simulation ---
print(f"Simulating FMU: {fmu_path}...")
try:
    # 1. Read Model Description
    model_description = read_model_description(fmu_path)

    # Prepare structured numpy array for input
    input_dtype = [('time', np.float64), ('u', np.float64, (13,))]
    input_data = np.zeros(2, dtype=input_dtype)
    input_data[0] = (0.0, u_val)
    input_data[1] = (t_stop, u_val)

    # We use CVode for better accuracy
    results = simulate_fmu(fmu_path, 
                          stop_time=t_stop, 
                          step_size=dt,
                          input=input_data,
                          output=['y'],
                          solver='CVode')
    
    time_list = results['time']
    # y is [Fx, Fy, Fz, Mx, My, Mz] in Body Frame
    # Since AeroPanels default is Z-UP body frame, Lift is +Fz
    Fz = results['y'][:, 2] 
    
    rho = 1.225
    S = 100.0
    CL_list = Fz / (0.5 * rho * V**2 * S)
    
    CL_steady_theory = 2 * np.pi * alpha_rad
    CL_normalized = CL_list / CL_steady_theory

    print(f"Final CL (Normalized): {CL_normalized[-1]:.4f}")
    print(f"Error vs Theory: {abs(CL_normalized[-1] - 1.0) * 100:.2f}%")

    # --- 3. Plot Results ---
    plt.figure(figsize=(10, 6))
    plt.plot(time_list, CL_normalized, label='FMU Export (FMI 3.0)')
    plt.axhline(1.0, color='r', linestyle='--', label='Steady Theory (2*pi*alpha)')
    plt.xlabel('Time [s]')
    plt.ylabel('CL / CL_steady')
    plt.title('Wagner Problem Verification (Normalized)')
    plt.grid(True)
    plt.ylim([0, 1.2])
    plt.legend()
    plt.savefig('examples/fmu_export/wagner_results.png')
    print("Simulation successful! Results saved to examples/fmu_export/wagner_results.png")

except Exception as e:
    print(f"Simulation failed: {e}")
