module FMI_XML

using XML
using AeroPanels

export generate_xml_node

function fmi2_real_variable(name, vr, causality; variability="continuous", initial=nothing, start=nothing, derivative=nothing)
    h = XML.h
    attrs = Dict{Symbol, Any}(
        :name => name,
        :valueReference => string(vr),
        :causality => causality,
        :variability => variability
    )
    if !isnothing(initial)
        attrs[:initial] = initial
    end
    
    var_node = h.ScalarVariable(; attrs...)
    
    real_attrs = Dict{Symbol, Any}()
    if !isnothing(start)
        real_attrs[:start] = string(start)
    end
    if !isnothing(derivative)
        real_attrs[:derivative] = string(derivative)
    end
    
    push!(var_node, h.Real(; real_attrs...))
    return var_node
end

"""
    generate_xml_node(model::AbstractFMUModel; modelName="...", version="1.0.0")

Generates the complete FMI 2.0 `modelDescription.xml` schema directly from `model` and its `@fmu_model` layout.
"""
function generate_xml_node(model::AbstractFMUModel; modelName::String=string(typeof(model).name.name), version::String="1.0.0")
    h = XML.h
    layout = GetFMULayout(model)
    arr = CreateFMUArray(model)
    
    root = h.fmiModelDescription(
        fmiVersion="2.0",
        modelName=modelName,
        guid="{$(bytes2hex(rand(UInt8, 16)))}",
        generationTool="AeroPanels.jl",
        version=version,
        variableNamingConvention="structured",
        numberOfEventIndicators="0"
    )
    
    push!(root, h.ModelExchange(
        modelIdentifier=modelName,
        canGetAndSetFMUstate="true",
        canSerializeFMUstate="false",
        providesDirectionalDerivative="true"
    ))
    
    push!(root, h.CoSimulation(
        modelIdentifier=modelName,
        canHandleVariableCommunicationStepSize="true",
        canGetAndSetFMUstate="true",
        canSerializeFMUstate="false",
        providesDirectionalDerivative="true"
    ))
    
    vars = h.ModelVariables()
    output_indices = Int[]
    derivative_indices = Int[]
    
    vr = 1
    
    # 1. States
    out_states = hasproperty(layout, :OutputStates) ? layout.OutputStates : ()
    for fname in layout.States
        sub_arr = getproperty(arr, fname)
        sz = length(sub_arr)
        is_out = fname in out_states
        causality_str = is_out ? "output" : "local"
        for i in 1:sz
            var_name = sz == 1 ? string(fname) : "$(fname)[$i]"
            push!(vars, fmi2_real_variable(var_name, vr, causality_str; initial="exact", start=0.0))
            if is_out
                push!(output_indices, vr)
            end
            vr += 1
        end
    end
    
    # 2. State Derivatives
    out_derivs = hasproperty(layout, :OutputDerivatives) ? layout.OutputDerivatives : ()
    for fname in layout.Derivatives
        sub_arr = getproperty(arr, fname)
        sz = length(sub_arr)
        is_out = fname in out_derivs
        causality_str = is_out ? "output" : "local"
        for i in 1:sz
            var_name = sz == 1 ? string(fname) : "$(fname)[$i]"
            state_vr = i # Maps 1-to-1 with corresponding state VR
            push!(vars, fmi2_real_variable(var_name, vr, causality_str; derivative=state_vr, initial="calculated"))
            push!(derivative_indices, vr)
            if is_out
                push!(output_indices, vr)
            end
            vr += 1
        end
    end
    
    # 3. Inputs
    for fname in layout.Inputs
        sub_arr = getproperty(arr, fname)
        sz = length(sub_arr)
        start_val = fname == :rho ? 1.225 : 0.0
        for i in 1:sz
            var_name = sz == 1 ? string(fname) : "$(fname)[$i]"
            push!(vars, fmi2_real_variable(var_name, vr, "input"; start=start_val))
            vr += 1
        end
    end
    
    # 4. Outputs
    mp_list = hasproperty(model, :monitorPoints) ? model.monitorPoints : MonitorPoint[]
    for fname in layout.Outputs
        sub_arr = getproperty(arr, fname)
        sz = length(sub_arr)
        if (fname == :Fmp || fname == :monitorPointLoads) && length(mp_list) > 0
            for (k, mp) in enumerate(mp_list)
                for c_idx in 1:6
                    var_name = "$(mp.name)_load[$c_idx]"
                    push!(vars, fmi2_real_variable(var_name, vr, "output"; initial="calculated"))
                    push!(output_indices, vr)
                    vr += 1
                end
            end
        else
            for i in 1:sz
                var_name = sz == 1 ? string(fname) : "$(fname)[$i]"
                push!(vars, fmi2_real_variable(var_name, vr, "output"; initial="calculated"))
                push!(output_indices, vr)
                vr += 1
            end
        end
    end
    
    push!(root, vars)
    
    # ModelStructure with optional Sparsity Dependencies
    sparsity_info = try
        GetFMISparsity(model)
    catch
        nothing
    end
    
    struct_node = h.ModelStructure()
    
    outputs_node = h.Outputs()
    if sparsity_info !== nothing
        J = sparsity_info.jacobian
        nStates = NumberOfStates(model)
        for (i, idx) in enumerate(output_indices)
            row_idx = nStates + i
            if row_idx <= size(J, 1)
                deps = findall(J[row_idx, :])
                if !isempty(deps)
                    deps_str = join(deps, " ")
                    kinds_str = join(["dependent" for _ in 1:length(deps)], " ")
                    push!(outputs_node, h.Unknown(index=string(idx), dependencies=deps_str, dependenciesKind=kinds_str))
                else
                    push!(outputs_node, h.Unknown(index=string(idx), dependencies="", dependenciesKind=""))
                end
            else
                push!(outputs_node, h.Unknown(index=string(idx)))
            end
        end
    else
        for idx in output_indices
            push!(outputs_node, h.Unknown(index=string(idx)))
        end
    end
    push!(struct_node, outputs_node)
    
    if length(derivative_indices) > 0
        derivatives_node = h.Derivatives()
        if sparsity_info !== nothing
            J = sparsity_info.jacobian
            nStates = NumberOfStates(model)
            for i in 1:length(derivative_indices)
                idx = derivative_indices[i]
                deps = findall(J[i, :])
                if !isempty(deps)
                    deps_str = join(deps, " ")
                    kinds = [d <= nStates ? "fixed" : "dependent" for d in deps]
                    kinds_str = join(kinds, " ")
                    push!(derivatives_node, h.Unknown(index=string(idx), dependencies=deps_str, dependenciesKind=kinds_str))
                else
                    push!(derivatives_node, h.Unknown(index=string(idx), dependencies="", dependenciesKind=""))
                end
            end
        else
            for idx in derivative_indices
                push!(derivatives_node, h.Unknown(index=string(idx)))
            end
        end
        push!(struct_node, derivatives_node)
    end
    
    initial_node = h.InitialUnknowns()
    for idx in output_indices
        push!(initial_node, h.Unknown(index=string(idx)))
    end
    for idx in derivative_indices
        push!(initial_node, h.Unknown(index=string(idx)))
    end
    push!(struct_node, initial_node)
    
    push!(root, struct_node)
    
    return root
end

end
