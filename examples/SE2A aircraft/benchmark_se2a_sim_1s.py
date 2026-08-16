import os
import time
import fmpy
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
import numpy as np

fmu_path = 'examples/SE2A aircraft/SE2AAircraftFMU.fmu'

if not os.path.exists(fmu_path):
    print(f"Error: {fmu_path} not found. Please run the Julia build script first.")
    exit(1)

print("Starting SE2A Aircraft FMU 1-Second Co-Simulation Benchmark...")
print(f"FMU Path: {fmu_path}")

overall_start = time.perf_counter()

# --- 1. FMU Loading and Instantiation ---
init_start = time.perf_counter()
try:
    unzipdir = extract(fmu_path)
    model_description = read_model_description(fmu_path)
    model_id = model_description.coSimulation.modelIdentifier

    def get_vrs(prefix, size):
        return [next(v.valueReference for v in model_description.modelVariables if v.name == f"{prefix}[{i}]") for i in range(1, size + 1)]

    vr_vb = get_vrs("vb", 3)
    vr_wb = get_vrs("wb", 3)
    vr_ab = get_vrs("ab", 3)
    vr_dwb = get_vrs("dwb", 3)
    
    vr_rho = [next(v.valueReference for v in model_description.modelVariables if v.name == "rho")]
    
    vr_PortElevator_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "u_cs.PortElevator")]
    vr_Rudder_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "u_cs.Rudder")]
    
    vr_rs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("rs[")]
    vr_vs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("vs[")]
    vr_as = [v.valueReference for v in model_description.modelVariables if v.name.startswith("as[")]

    vr_Fb = get_vrs("Fb", 3)
    vr_Mb = get_vrs("Mb", 3)
    
    vr_StbdWingRoot_load = get_vrs("StbdWingRoot_load", 6)
    vr_PortWingRoot_load = get_vrs("PortWingRoot_load", 6)

    fmu = FMU2Slave(
        guid=model_description.guid,
        unzipDirectory=unzipdir,
        instanceName='instance1',
        modelIdentifier=model_id
    )

    fmu.instantiate()

    t_stop = 1.0  # 1 second simulation
    dt = 0.01     # 100 Hz communication interval (100 steps total)
    vb = [100.0, 0.0, 10.0]

    fmu.setupExperiment(startTime=0.0)
    fmu.enterInitializationMode()
    fmu.setReal(vr_vb, vb)
    fmu.setReal(vr_wb, [0.0]*3)
    fmu.setReal(vr_ab, [0.0]*3)
    fmu.setReal(vr_dwb, [0.0]*3)
    fmu.setReal(vr_rho, [1.225])
    
    fmu.setReal(vr_PortElevator_delta, [-0.05])
    fmu.setReal(vr_Rudder_delta, [0.1])
    
    if vr_rs:
        fmu.setReal(vr_rs, [0.0]*len(vr_rs))
    if vr_vs:
        fmu.setReal(vr_vs, [0.0]*len(vr_vs))
    if vr_as:
        fmu.setReal(vr_as, [0.0]*len(vr_as))
        
    fmu.exitInitializationMode()
    init_end = time.perf_counter()
    
    # --- 2. Simulation Loop ---
    print("\nRunning co-simulation loop...")
    sim_start = time.perf_counter()
    
    current_time = 0.0
    steps = 0
    while current_time <= t_stop + 1e-9:
        # Retrieve outputs at each step
        Fb = fmu.getReal(vr_Fb)
        Mb = fmu.getReal(vr_Mb)
        load_stbd = fmu.getReal(vr_StbdWingRoot_load)
        load_port = fmu.getReal(vr_PortWingRoot_load)
        
        if current_time >= t_stop - 1e-9:
            break
            
        fmu.doStep(currentCommunicationPoint=current_time, communicationStepSize=dt)
        current_time += dt
        steps += 1
        
    sim_end = time.perf_counter()
    
    fmu.terminate()
    fmu.freeInstance()
    
    overall_end = time.perf_counter()
    
    init_time = init_end - init_start
    sim_time = sim_end - sim_start
    total_time = overall_end - overall_start
    
    print("\n================ CO-SIMULATION BENCHMARK RESULTS ================")
    print(f"Simulation Duration:     {t_stop:.1f} seconds")
    print(f"Communication Step (dt): {dt:.3f} seconds ({steps} steps total)")
    print(f"Loading & Init Time:     {init_time:.4f} seconds")
    print(f"Pure Simulation Time:    {sim_time:.4f} seconds")
    print(f"Total Execution Time:    {total_time:.4f} seconds")
    print(f"Average Time per Step:   {(sim_time / steps) * 1000:.2f} milliseconds")
    print("==================================================================")
    
    print("\nLast computed state outputs at t=1.0s:")
    print(f"  Fb (Forces): {Fb}")
    print(f"  Mb (Moments): {Mb}")
    print(f"  StbdWingRoot loads: {load_stbd}")
    print(f"  PortWingRoot loads: {load_port}")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Co-Simulation Benchmark failed: {e}")
