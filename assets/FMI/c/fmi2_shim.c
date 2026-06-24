#include <windows.h>
#include "fmi2Functions.h"
#include <string.h>

/*
 * GLOBAL VARIABLES
 * g_hinst: Stores the instance handle of THIS shim DLL.
 * julia_dll: Stores the handle to the loaded Julia engine DLL.
 */
static HINSTANCE g_hinst = NULL;
static HMODULE julia_dll = NULL;

/*
 * DllMain is the entry point for Windows DLLs.
 * We use it to capture the instance handle (g_hinst) as soon as the DLL is loaded
 * by the importer.
 */
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved) {
    if (fdwReason == DLL_PROCESS_ATTACH) {
        g_hinst = hinstDLL;
        DisableThreadLibraryCalls(hinstDLL);
    }
    return TRUE;
}

/*
 * load_julia_engine() solves the Windows DLL search path problem.
 */
static void load_julia_engine() {
    if (julia_dll != NULL) return;

    char path[MAX_PATH];
    GetModuleFileNameA(g_hinst, path, MAX_PATH);

    char *last_slash = strrchr(path, '\\');
    if (last_slash) *last_slash = '\0';

    SetDllDirectoryA(path);
    julia_dll = LoadLibraryA("fmu_engine.dll");
}

/*
 * MACRO: FMI2_FWD
 * Generates the boilerplate wrapper for any FMI 2.0 function.
 */
#define FMI2_FWD(ret, name, sig, args, err_ret) \
    FMI2_Export ret name sig { \
        load_julia_engine(); \
        if (!julia_dll) return err_ret; \
        typedef ret (*ptr_t) sig; \
        ptr_t func = (ptr_t)GetProcAddress(julia_dll, #name); \
        if (!func) return err_ret; \
        return func args; \
    }

#define FMI2_STATUS(name, sig, args) FMI2_FWD(fmi2Status, name, sig, args, fmi2Error)

/* =========================================================================
 * FMI 2.0 FUNCTION FORWARDERS
 * ========================================================================= */

FMI2_FWD(const char*, fmi2GetTypesPlatform,
         (void),
         (),
         "default")

FMI2_FWD(const char*, fmi2GetVersion,
         (void),
         (),
         "2.0")

FMI2_STATUS(fmi2SetDebugLogging,
            (fmi2Component c, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]),
            (c, loggingOn, nCategories, categories))

FMI2_FWD(fmi2Component, fmi2Instantiate,
         (fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGUID, fmi2String fmuResourceLocation, const fmi2CallbackFunctions* functions, fmi2Boolean visible, fmi2Boolean loggingOn),
         (instanceName, fmuType, fmuGUID, fmuResourceLocation, functions, visible, loggingOn),
         NULL)

FMI2_FWD(void, fmi2FreeInstance,
         (fmi2Component c),
         (c), )

FMI2_STATUS(fmi2SetupExperiment,
            (fmi2Component c, fmi2Boolean toleranceDefined, fmi2Real tolerance, fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime),
            (c, toleranceDefined, tolerance, startTime, stopTimeDefined, stopTime))

FMI2_STATUS(fmi2EnterInitializationMode,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2ExitInitializationMode,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2Terminate,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2Reset,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2GetReal,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2GetInteger,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2GetBoolean,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2GetString,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2SetReal,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2SetInteger,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2SetBoolean,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2SetString,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[]),
            (c, vr, nvr, value))

FMI2_STATUS(fmi2GetFMUstate,
            (fmi2Component c, fmi2FMUstate* FMUstate),
            (c, FMUstate))

FMI2_STATUS(fmi2SetFMUstate,
            (fmi2Component c, fmi2FMUstate FMUstate),
            (c, FMUstate))

FMI2_STATUS(fmi2FreeFMUstate,
            (fmi2Component c, fmi2FMUstate* FMUstate),
            (c, FMUstate))

FMI2_STATUS(fmi2SerializedFMUstateSize,
            (fmi2Component c, fmi2FMUstate FMUstate, size_t* size),
            (c, FMUstate, size))

FMI2_STATUS(fmi2SerializeFMUstate,
            (fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size),
            (c, FMUstate, serializedState, size))

FMI2_STATUS(fmi2DeSerializeFMUstate,
            (fmi2Component c, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate),
            (c, serializedState, size, FMUstate))

FMI2_STATUS(fmi2GetDirectionalDerivative,
            (fmi2Component c, const fmi2ValueReference vUnknown_ref[], size_t nUnknown, const fmi2ValueReference vKnown_ref[], size_t nKnown, const fmi2Real dvKnown[], fmi2Real dvUnknown[]),
            (c, vUnknown_ref, nUnknown, vKnown_ref, nKnown, dvKnown, dvUnknown))

/* Model Exchange Functions */
FMI2_STATUS(fmi2EnterEventMode,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2NewDiscreteStates,
            (fmi2Component c, fmi2EventInfo* eventInfo),
            (c, eventInfo))

FMI2_STATUS(fmi2EnterContinuousTimeMode,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2CompletedIntegratorStep,
            (fmi2Component c, fmi2Boolean noSetFMUStatePriorToCurrentPoint, fmi2Boolean* enterEventMode, fmi2Boolean* terminateSimulation),
            (c, noSetFMUStatePriorToCurrentPoint, enterEventMode, terminateSimulation))

FMI2_STATUS(fmi2SetTime,
            (fmi2Component c, fmi2Real time),
            (c, time))

FMI2_STATUS(fmi2SetContinuousStates,
            (fmi2Component c, const fmi2Real x[], size_t nx),
            (c, x, nx))

FMI2_STATUS(fmi2GetDerivatives,
            (fmi2Component c, fmi2Real dx[], size_t nx),
            (c, dx, nx))

FMI2_STATUS(fmi2GetEventIndicators,
            (fmi2Component c, fmi2Real eventIndicators[], size_t ni),
            (c, eventIndicators, ni))

FMI2_STATUS(fmi2GetContinuousStates,
            (fmi2Component c, fmi2Real x[], size_t nx),
            (c, x, nx))

FMI2_STATUS(fmi2GetNominalsOfContinuousStates,
            (fmi2Component c, fmi2Real x_nominal[], size_t nx),
            (c, x_nominal, nx))

/* Co-Simulation Functions */
FMI2_STATUS(fmi2SetRealInputDerivatives,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], const fmi2Real value[]),
            (c, vr, nvr, order, value))

FMI2_STATUS(fmi2GetRealOutputDerivatives,
            (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], fmi2Real value[]),
            (c, vr, nvr, order, value))

FMI2_STATUS(fmi2DoStep,
            (fmi2Component c, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize, fmi2Boolean noSetFMUStatePriorToCurrentPoint),
            (c, currentCommunicationPoint, communicationStepSize, noSetFMUStatePriorToCurrentPoint))

FMI2_STATUS(fmi2CancelStep,
            (fmi2Component c),
            (c))

FMI2_STATUS(fmi2GetStatus,
            (fmi2Component c, const fmi2StatusKind s, fmi2Status* value),
            (c, s, value))

FMI2_STATUS(fmi2GetRealStatus,
            (fmi2Component c, const fmi2StatusKind s, fmi2Real* value),
            (c, s, value))

FMI2_STATUS(fmi2GetIntegerStatus,
            (fmi2Component c, const fmi2StatusKind s, fmi2Integer* value),
            (c, s, value))

FMI2_STATUS(fmi2GetBooleanStatus,
            (fmi2Component c, const fmi2StatusKind s, fmi2Boolean* value),
            (c, s, value))

FMI2_STATUS(fmi2GetStringStatus,
            (fmi2Component c, const fmi2StatusKind s, fmi2String* value),
            (c, s, value))
