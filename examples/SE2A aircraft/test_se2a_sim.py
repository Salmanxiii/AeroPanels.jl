import os
import fmpy
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
import numpy as np

fmu_path = 'examples/SE2A aircraft/SE2AAircraftFMU.fmu'

if not os.path.exists(fmu_path):
    print(f"Error: {fmu_path} not found. Please run the Julia build script first.")
    exit(1)

print(f"Simulating SE2A Aircraft FMU: {fmu_path}...")
try:
    # 1. Read Model Description and Extract FMU
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
    vr_PortElevator_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "PortElevator_delta")]
    vr_Rudder_delta = [next(v.valueReference for v in model_description.modelVariables if v.name == "Rudder_delta")]
    
    # Optional structural nodes
    vr_rs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("rs[")]
    vr_vs = [v.valueReference for v in model_description.modelVariables if v.name.startswith("vs[")]
    vr_as = [v.valueReference for v in model_description.modelVariables if v.name.startswith("as[")]

    # Outputs
    vr_Fb = get_vrs("Fb", 3)
    vr_Mb = get_vrs("Mb", 3)
    
    vr_StbdWingRoot_load = get_vrs("StbdWingRoot_load", 6)
    vr_PortWingRoot_load = get_vrs("PortWingRoot_load", 6)

    # Instantiate FMU2Slave
    fmu = FMU2Slave(
        guid=model_description.guid,
        unzipDirectory=unzipdir,
        instanceName='instance1',
        modelIdentifier=model_id
    )

    fmu.instantiate()

    # Define simulation values
    t_stop = 1.0
    dt = 0.1
    V = 30.0  # 30 m/s forward speed
    vb = [V, 0.0, 0.0]

    # Enter initialization mode and set inputs
    fmu.setupExperiment(startTime=0.0)
    fmu.enterInitializationMode()
    fmu.setReal(vr_vb, vb)
    fmu.setReal(vr_wb, [0.0]*3)
    fmu.setReal(vr_ab, [0.0]*3)
    fmu.setReal(vr_dwb, [0.0]*3)
    fmu.setReal(vr_rho, [1.225])
    
    # Deflect control surfaces
    fmu.setReal(vr_PortElevator_delta, [-0.05])  # -0.05 rad elevator
    fmu.setReal(vr_Rudder_delta, [0.1])          # 0.1 rad rudder
    
    if vr_rs:
        fmu.setReal(vr_rs, [0.0]*len(vr_rs))
    if vr_vs:
        fmu.setReal(vr_vs, [0.0]*len(vr_vs))
    if vr_as:
        fmu.setReal(vr_as, [0.0]*len(vr_as))
        
    fmu.exitInitializationMode()

    # Simulation loop
    current_time = 0.0
    while current_time <= t_stop + 1e-9:
        Fb = fmu.getReal(vr_Fb)
        Mb = fmu.getReal(vr_Mb)
        load_stbd = fmu.getReal(vr_StbdWingRoot_load)
        load_port = fmu.getReal(vr_PortWingRoot_load)
        
        print(f"t={current_time:.1f}s:")
        print(f"  Fb (Forces): {Fb}")
        print(f"  Mb (Moments): {Mb}")
        print(f"  StbdWingRoot loads (Force & Moment): {load_stbd}")
        print(f"  PortWingRoot loads (Force & Moment): {load_port}")
        
        if current_time >= t_stop - 1e-9:
            break
            
        fmu.doStep(currentCommunicationPoint=current_time, communicationStepSize=dt)
        current_time += dt

    fmu.terminate()
    fmu.freeInstance()
    print("\nSimulation and verification completed successfully!")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Simulation failed: {e}")
