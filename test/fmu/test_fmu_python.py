import os
import time
from fmpy import read_model_description, instantiate_fmu, extract

def test_prebuilt_fmu():
    fmu_path = os.path.abspath(r"scratch\prebuilt_fmu_test\SE2APrebuiltSteadyFMU.fmu")
    print(f"Testing prebuilt FMU with FMPy: {fmu_path}")
    assert os.path.exists(fmu_path), f"FMU binary not found at {fmu_path}"
    
    model_description = read_model_description(fmu_path)
    unzipdir = extract(fmu_path)
    
    t0 = time.time()
    fmu = instantiate_fmu(unzipdir, model_description, fmi_type='CoSimulation')
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()
    t1 = time.time()
    print(f"FMU Instantation and Initialization time: {(t1 - t0)*1000:.2f} ms")
    
    fmu.doStep(currentCommunicationPoint=0.0, communicationStepSize=0.01)
    
    vrs = {v.name: v.valueReference for v in model_description.modelVariables}
    print("Variables in FMU:")
    for name in sorted(vrs.keys()):
        if "lift" in name.lower() or "drag" in name.lower() or "moment" in name.lower() or "cl" in name.lower():
            val = fmu.getReal([vrs[name]])[0]
            print(f"  {name} = {val:.6f}")
    
    fmu.terminate()
    fmu.freeInstance()
    print("SUCCESS: Python FMPy simulation of prebuilt FMU finished cleanly!")

if __name__ == "__main__":
    test_prebuilt_fmu()
