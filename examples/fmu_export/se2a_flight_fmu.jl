using AeroPanels

println("Building 6DOF Flight Dynamics FMU for SE2A Aircraft...")
flight_model = SE2AUnsteadyFlight(; nWake=40, wakeLength=10.0);

fmu_path = BuildFMU(flight_model; fmu_name="SE2AFlightFMU", output_dir=@__DIR__)
println("Successfully created 6DOF Flight Dynamics FMU at: ", fmu_path)
