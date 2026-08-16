import os
import time
import fmpy
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Model
import numpy as np
from scipy.integrate import solve_ivp
from ctypes import POINTER, c_double

fmu_path = 'examples/SE2A aircraft/SE2AAircraftFMU.fmu'

if not os.path.exists(fmu_path):
    print(f"Error: {fmu_path} not found. Please run the Julia build script first.")
    exit(1)

print("Starting SE2A Aircraft FMU Model Exchange Benchmark...")
print(f"FMU Path: {fmu_path}")

overall_start = time.perf_counter()

# --- 1. FMU Loading and Instantiation ---
init_start = time.perf_counter()
try:
    unzipdir = extract(fmu_path)
    model_description = read_model_description(fmu_path)
    model_id = model_description.coSimulation.modelIdentifier if model_description.coSimulation is not None else model_description.modelExchange.modelIdentifier

    # Lookup value references dynamically
    def get_vrs(prefix, size):
        return [next(v.valueReference for v in model_description.modelVariables if v.name == f"{prefix}[{i}]") for i in range(1, size + 1)]

    vr_vb = get_vrs("vb", 3)
    vr_wb = get_vrs("wb", 3)
    vr_ab = get_vrs("ab", 3)
    vr_dwb = get_vrs("dwb", 3)
    
    vr_rho = [next(v.valueReference for v in model_description.modelVariables if v.name == "rho")]
    
    # Control surfaces
    vr_PortElevator_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "u_cs.PortElevator")]
    vr_Rudder_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "u_cs.Rudder")]
    
    # Optional structural nodes
    vr_rs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("rs[")]
    vr_vs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("vs[")]
    vr_as = [v.valueReference for v in model_description.modelVariables if v.name.startswith("as[")]

    # Outputs
    vr_Fb = get_vrs("Fb", 3)
    vr_Mb = get_vrs("Mb", 3)
    
    vr_StbdWingRoot_load = get_vrs("StbdWingRoot_load", 6)
    vr_PortWingRoot_load = get_vrs("PortWingRoot_load", 6)

    # Instantiate FMU2Model (Model Exchange)
    fmu = FMU2Model(
        guid=model_description.guid,
        unzipDirectory=unzipdir,
        instanceName='instance1',
        modelIdentifier=model_id
    )

    fmu.instantiate()

    # Define simulation values
    t_stop = 1.0  # 1 second simulation
    vb = [100.0, 0.0, 10.0]

    # Enter initialization mode and set inputs
    fmu.setupExperiment(startTime=0.0)
    fmu.enterInitializationMode()
    fmu.setReal(vr_vb, vb)
    fmu.setReal(vr_wb, [0.0]*3)
    fmu.setReal(vr_ab, [0.0]*3)
    fmu.setReal(vr_dwb, [0.0]*3)
    fmu.setReal(vr_rho, [1.225])
    
    # Deflect control surfaces
    fmu.setReal(vr_PortElevator_delta, [-0.05])
    fmu.setReal(vr_Rudder_delta, [0.1])
    
    if vr_rs:
        fmu.setReal(vr_rs, [0.0]*len(vr_rs))
    if vr_vs:
        fmu.setReal(vr_vs, [0.0]*len(vr_vs))
    if vr_as:
        fmu.setReal(vr_as, [0.0]*len(vr_as))
        
    fmu.exitInitializationMode()
    fmu.enterContinuousTimeMode()
    init_end = time.perf_counter()
    
    # --- 2. Model Exchange ODE Definition ---
    # Number of states
    nx = model_description.numberOfContinuousStates
    print(f"Continuous States (nx): {nx}")

    # Pre-allocate derivative array and its ctypes pointer to avoid memory allocations in loop
    dx = np.zeros(nx, dtype=np.float64)
    pdx = dx.ctypes.data_as(POINTER(c_double))

    # Helper function to compute derivatives dy/dt
    evaluations = 0
    def ode_func(t, y):
        global evaluations
        evaluations += 1
        
        # Set time in the FMU
        fmu.setTime(t)
        
        # Convert y to ctypes pointer and set continuous states in the FMU
        # y is a 1D numpy array passed by solve_ivp
        py = y.ctypes.data_as(POINTER(c_double))
        fmu.setContinuousStates(py, nx)
        
        # Retrieve derivatives into pre-allocated dx
        fmu.getDerivatives(pdx, nx)
        
        return dx

    # --- 3. Run Variable-Step Simulation ---
    print("\nRunning Model Exchange simulation using SciPy solve_ivp (RK45 - Dormand-Prince)...")
    sim_start = time.perf_counter()
    
    # Initial state vector (all zeros)
    y0 = np.zeros(nx, dtype=np.float64)
    
    # Solve ODE using variable step Runge-Kutta 45 (Dormand-Prince)
    sol = solve_ivp(
        ode_func,
        (0.0, t_stop),
        y0,
        method='RK45',
        rtol=1e-3,
        atol=1e-6
    )
    
    sim_end = time.perf_counter()
    
    # --- 4. Retrieve Outputs at Final Step ---
    # Set final state to get outputs
    fmu.setTime(sol.t[-1])
    
    y_final = np.ascontiguousarray(sol.y[:, -1], dtype=np.float64)
    py_final = y_final.ctypes.data_as(POINTER(c_double))
    fmu.setContinuousStates(py_final, nx)
    
    Fb = fmu.getReal(vr_Fb)
    Mb = fmu.getReal(vr_Mb)
    load_stbd = fmu.getReal(vr_StbdWingRoot_load)
    load_port = fmu.getReal(vr_PortWingRoot_load)
    
    # Terminate and free FMU
    fmu.terminate()
    fmu.freeInstance()
    
    overall_end = time.perf_counter()
    
    # --- 5. Report Benchmark Results ---
    init_time = init_end - init_start
    sim_time = sim_end - sim_start
    total_time = overall_end - overall_start
    
    print("\n================ BENCHMARK RESULTS (MODEL EXCHANGE) ================")
    print(f"Simulation Duration:     {t_stop:.1f} seconds")
    print(f"Solver Method:           RK45 (Dormand-Prince, variable-step)")
    print(f"Successful Steps:        {sol.t.size} steps")
    print(f"Function Evaluations:    {evaluations} calls")
    print(f"Loading & Init Time:     {init_time:.4f} seconds")
    print(f"Pure Simulation Time:    {sim_time:.4f} seconds")
    print(f"Total Execution Time:    {total_time:.4f} seconds")
    print(f"Average Time per Call:   {(sim_time / evaluations) * 1000:.2f} milliseconds")
    print("=====================================================================")
    
    print("\nLast computed state outputs at t=10.0s:")
    print(f"  Fb (Forces): {Fb}")
    print(f"  Mb (Moments): {Mb}")
    print(f"  StbdWingRoot loads: {load_stbd}")
    print(f"  PortWingRoot loads: {load_port}")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"ME Benchmark failed: {e}")
