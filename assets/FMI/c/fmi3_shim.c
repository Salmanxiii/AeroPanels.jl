#include <windows.h>
#include <stdio.h>
#include <string.h>
#include "fmi3Functions.h"

/* 
 * GLOBAL VARIABLE
 * julia_dll: Stores the handle to the loaded Julia engine DLL.
 */
static HMODULE julia_dll = NULL;

/*
 * load_julia_engine() solves the Windows DLL search path problem.
 * It uses GetModuleHandleExA to find its own location robustly.
 */
static void load_julia_engine() {
    if (julia_dll != NULL) return;
    
    HMODULE hShim = NULL;
    // GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS = 0x00000004
    // GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT = 0x00000002
    if (GetModuleHandleExA(0x00000004 | 0x00000002, (LPCSTR)&load_julia_engine, &hShim)) {
        char path[MAX_PATH];
        if (GetModuleFileNameA(hShim, path, MAX_PATH)) {
            char *last_slash = strrchr(path, '\\');
            if (last_slash) {
                *last_slash = '\0';
                SetDllDirectoryA(path);
                julia_dll = LoadLibraryA("fmu_engine.dll");
            }
        }
    }
}

/*
 * MACRO: FMI3_FWD
 * Generates the boilerplate wrapper for any FMI 3.0 function.
 * This makes the shim a Universal Forwarder for ME, CS, and SE.
 */
#define FMI3_FWD(ret, name, sig, args, err_ret) \
    FMI3_Export ret name sig { \
        load_julia_engine(); \
        if (!julia_dll) return err_ret; \
        typedef ret (*ptr_t) sig; \
        ptr_t func = (ptr_t)GetProcAddress(julia_dll, #name); \
        if (!func) return err_ret; \
        return func args; \
    }

#define FMI3_STATUS(name, sig, args) FMI3_FWD(fmi3Status, name, sig, args, fmi3Error)

/* =========================================================================
 * FMI 3.0 FUNCTION FORWARDERS (Universal Switchboard)
 * ========================================================================= */

FMI3_FWD(const char*, fmi3GetVersion, (void), (), NULL)

FMI3_STATUS(fmi3SetDebugLogging, 
            (fmi3Instance instance, fmi3Boolean loggingOn, size_t nCategories, const fmi3String categories[]), 
            (instance, loggingOn, nCategories, categories))

FMI3_FWD(fmi3Instance, fmi3InstantiateModelExchange, 
         (fmi3String instanceName, fmi3String instantiationToken, fmi3String resourcePath, 
          fmi3Boolean visible, fmi3Boolean loggingOn, fmi3InstanceEnvironment instanceEnvironment, 
          fmi3LogMessageCallback logMessage), 
         (instanceName, instantiationToken, resourcePath, visible, loggingOn, instanceEnvironment, logMessage), 
         NULL)

FMI3_FWD(fmi3Instance, fmi3InstantiateCoSimulation, 
         (fmi3String instanceName, fmi3String instantiationToken, fmi3String resourcePath, 
          fmi3Boolean visible, fmi3Boolean loggingOn, fmi3Boolean eventModeUsed, 
          fmi3Boolean earlyReturnAllowed, const fmi3ValueReference requiredIntermediateVariables[], 
          size_t nRequiredIntermediateVariables, fmi3InstanceEnvironment instanceEnvironment, 
          fmi3LogMessageCallback logMessage, fmi3IntermediateUpdateCallback intermediateUpdate), 
         (instanceName, instantiationToken, resourcePath, visible, loggingOn, eventModeUsed, earlyReturnAllowed, requiredIntermediateVariables, nRequiredIntermediateVariables, instanceEnvironment, logMessage, intermediateUpdate), 
         NULL)

FMI3_FWD(fmi3Instance, fmi3InstantiateScheduledExecution, 
         (fmi3String instanceName, fmi3String instantiationToken, fmi3String resourcePath, 
          fmi3Boolean visible, fmi3Boolean loggingOn, fmi3InstanceEnvironment instanceEnvironment, 
          fmi3LogMessageCallback logMessage, fmi3ClockUpdateCallback clockUpdate, 
          fmi3LockPreemptionCallback lockPreemption, fmi3UnlockPreemptionCallback unlockPreemption), 
         (instanceName, instantiationToken, resourcePath, visible, loggingOn, instanceEnvironment, logMessage, clockUpdate, lockPreemption, unlockPreemption), 
         NULL)

FMI3_FWD(void, fmi3FreeInstance, (fmi3Instance instance), (instance), )

FMI3_STATUS(fmi3EnterInitializationMode, 
            (fmi3Instance instance, fmi3Boolean toleranceDefined, fmi3Float64 tolerance, 
             fmi3Float64 startTime, fmi3Boolean stopTimeDefined, fmi3Float64 stopTime), 
            (instance, toleranceDefined, tolerance, startTime, stopTimeDefined, stopTime))

FMI3_STATUS(fmi3ExitInitializationMode, (fmi3Instance instance), (instance))
FMI3_STATUS(fmi3EnterEventMode, (fmi3Instance instance), (instance))
FMI3_STATUS(fmi3Terminate, (fmi3Instance instance), (instance))
FMI3_STATUS(fmi3Reset, (fmi3Instance instance), (instance))

/* Data I/O */
FMI3_STATUS(fmi3GetFloat64, (fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr, fmi3Float64 values[], size_t nValues), (instance, vr, nvr, values, nValues))
FMI3_STATUS(fmi3SetFloat64, (fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr, const fmi3Float64 values[], size_t nValues), (instance, vr, nvr, values, nValues))

/* Continuous Time (Model Exchange) */
FMI3_STATUS(fmi3EnterContinuousTimeMode, (fmi3Instance instance), (instance))
FMI3_STATUS(fmi3SetTime, (fmi3Instance instance, fmi3Float64 time), (instance, time))
FMI3_STATUS(fmi3SetContinuousStates, (fmi3Instance instance, const fmi3Float64 continuousStates[], size_t nContinuousStates), (instance, continuousStates, nContinuousStates))
FMI3_STATUS(fmi3GetContinuousStateDerivatives, (fmi3Instance instance, fmi3Float64 derivatives[], size_t nContinuousStates), (instance, derivatives, nContinuousStates))
FMI3_STATUS(fmi3GetContinuousStates, (fmi3Instance instance, fmi3Float64 continuousStates[], size_t nContinuousStates), (instance, continuousStates, nContinuousStates))
FMI3_STATUS(fmi3GetNominalsOfContinuousStates, (fmi3Instance instance, fmi3Float64 nominals[], size_t nContinuousStates), (instance, nominals, nContinuousStates))
FMI3_STATUS(fmi3GetNumberOfContinuousStates, (fmi3Instance instance, size_t* nContinuousStates), (instance, nContinuousStates))
FMI3_STATUS(fmi3CompletedIntegratorStep, (fmi3Instance instance, fmi3Boolean noSetFMUStatePriorToCurrentPoint, fmi3Boolean* enterEventMode, fmi3Boolean* terminateSimulation), (instance, noSetFMUStatePriorToCurrentPoint, enterEventMode, terminateSimulation))

/* Discrete States & Events */
FMI3_STATUS(fmi3UpdateDiscreteStates, (fmi3Instance instance, fmi3Boolean* discreteStatesNeedUpdate, fmi3Boolean* terminateSimulation, fmi3Boolean* nominalsOfContinuousStatesChanged, fmi3Boolean* valuesOfContinuousStatesChanged, fmi3Boolean* nextEventTimeDefined, fmi3Float64* nextEventTime), (instance, discreteStatesNeedUpdate, terminateSimulation, nominalsOfContinuousStatesChanged, valuesOfContinuousStatesChanged, nextEventTimeDefined, nextEventTime))
FMI3_STATUS(fmi3GetNumberOfEventIndicators, (fmi3Instance instance, size_t* nEventIndicators), (instance, nEventIndicators))
FMI3_STATUS(fmi3GetEventIndicators, (fmi3Instance instance, fmi3Float64 eventIndicators[], size_t nEventIndicators), (instance, eventIndicators, nEventIndicators))

/* Co-Simulation */
FMI3_STATUS(fmi3EnterStepMode, (fmi3Instance instance), (instance))
FMI3_STATUS(fmi3DoStep, (fmi3Instance instance, fmi3Float64 currentCommunicationPoint, fmi3Float64 communicationStepSize, fmi3Boolean noSetFMUStatePriorToCurrentPoint, fmi3Boolean* eventHandlingNeeded, fmi3Boolean* terminateSimulation, fmi3Boolean* earlyReturn, fmi3Float64* lastSuccessfulTime), (instance, currentCommunicationPoint, communicationStepSize, noSetFMUStatePriorToCurrentPoint, eventHandlingNeeded, terminateSimulation, earlyReturn, lastSuccessfulTime))

/* Dummies for unimplemented types to keep standard-compliant */
FMI3_STATUS(fmi3GetFloat32, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Float32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetInt8,    (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Int8 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetUInt8,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt8 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetInt16,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Int16 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetUInt16,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt16 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetInt32,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Int32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetUInt32,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetInt64,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Int64 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetUInt64,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt64 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetBoolean, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Boolean v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetString,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3String v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3GetBinary,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, size_t vs[], fmi3Binary v[], size_t nv), (i, vr, nvr, vs, v, nv))
FMI3_STATUS(fmi3GetClock,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Clock v[]), (i, vr, nvr, v))
FMI3_STATUS(fmi3SetFloat32, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Float32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetInt8,    (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Int8 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetUInt8,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt8 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetInt16,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Int16 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetUInt16,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt16 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetInt32,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Int32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetUInt32,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt32 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetInt64,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Int64 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetUInt64,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt64 v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetBoolean, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Boolean v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetString,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3String v[], size_t nv), (i, vr, nvr, v, nv))
FMI3_STATUS(fmi3SetBinary,  (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const size_t vs[], const fmi3Binary v[], size_t nv), (i, vr, nvr, vs, v, nv))
FMI3_STATUS(fmi3SetClock,   (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Clock v[]), (i, vr, nvr, v))
FMI3_STATUS(fmi3GetFMUState, (fmi3Instance i, fmi3FMUState* s), (i, s))
FMI3_STATUS(fmi3SetFMUState, (fmi3Instance i, fmi3FMUState s), (i, s))
FMI3_STATUS(fmi3FreeFMUState, (fmi3Instance i, fmi3FMUState* s), (i, s))
FMI3_STATUS(fmi3SerializedFMUStateSize, (fmi3Instance i, fmi3FMUState s, size_t* sz), (i, s, sz))
FMI3_STATUS(fmi3SerializeFMUState, (fmi3Instance i, fmi3FMUState s, fmi3Byte st[], size_t sz), (i, s, st, sz))
FMI3_STATUS(fmi3DeserializeFMUState, (fmi3Instance i, const fmi3Byte st[], size_t sz, fmi3FMUState* s), (i, st, sz, s))
FMI3_STATUS(fmi3GetDirectionalDerivative, (fmi3Instance i, const fmi3ValueReference u[], size_t nu, const fmi3ValueReference k[], size_t nk, const fmi3Float64 sd[], size_t nsd, fmi3Float64 sens[], size_t nsens), (i, u, nu, k, nk, sd, nsd, sens, nsens))
FMI3_STATUS(fmi3GetAdjointDerivative, (fmi3Instance i, const fmi3ValueReference u[], size_t nu, const fmi3ValueReference k[], size_t nk, const fmi3Float64 sd[], size_t nsd, fmi3Float64 sens[], size_t nsens), (i, u, nu, k, nk, sd, nsd, sens, nsens))
FMI3_STATUS(fmi3EnterConfigurationMode, (fmi3Instance i), (i))
FMI3_STATUS(fmi3ExitConfigurationMode, (fmi3Instance i), (i))
FMI3_STATUS(fmi3GetIntervalDecimal, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Float64 intv[], fmi3IntervalQualifier q[]), (i, vr, nvr, intv, q))
FMI3_STATUS(fmi3GetIntervalFraction, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt64 ic[], fmi3UInt64 r[], fmi3IntervalQualifier q[]), (i, vr, nvr, ic, r, q))
FMI3_STATUS(fmi3GetShiftDecimal, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3Float64 s[]), (i, vr, nvr, s))
FMI3_STATUS(fmi3GetShiftFraction, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, fmi3UInt64 sc[], fmi3UInt64 r[]), (i, vr, nvr, sc, r))
FMI3_STATUS(fmi3SetIntervalDecimal, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Float64 intv[]), (i, vr, nvr, intv))
FMI3_STATUS(fmi3SetIntervalFraction, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt64 ic[], const fmi3UInt64 r[]), (i, vr, nvr, ic, r))
FMI3_STATUS(fmi3SetShiftDecimal, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Float64 s[]), (i, vr, nvr, s))
FMI3_STATUS(fmi3SetShiftFraction, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3UInt64 sc[], const fmi3UInt64 r[]), (i, vr, nvr, sc, r))
FMI3_STATUS(fmi3EvaluateDiscreteStates, (fmi3Instance i), (i))
FMI3_STATUS(fmi3GetOutputDerivatives, (fmi3Instance i, const fmi3ValueReference vr[], size_t nvr, const fmi3Int32 o[], fmi3Float64 v[], size_t nv), (i, vr, nvr, o, v, nv))
FMI3_STATUS(fmi3ActivateModelPartition, (fmi3Instance i, fmi3ValueReference clk, fmi3Float64 t), (i, clk, t))
FMI3_STATUS(fmi3GetNumberOfVariableDependencies, (fmi3Instance i, fmi3ValueReference vr, size_t* nd), (i, vr, nd))
FMI3_STATUS(fmi3GetVariableDependencies, (fmi3Instance i, fmi3ValueReference d, size_t edi[], fmi3ValueReference ids[], size_t eii[], fmi3DependencyKind dk[], size_t nd), (i, d, edi, ids, eii, dk, nd))
