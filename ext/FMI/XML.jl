module FMI_XML

using XML

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

function generate_xml_node(name, nx, nCtrl, csNames, nVert, nmp, mpNames, version)
    h = XML.h
    
    root = h.fmiModelDescription(
        fmiVersion="2.0",
        modelName=name,
        guid="{$(bytes2hex(rand(UInt8, 16)))}",
        generationTool="AeroPanels.jl",
        version=version,
        variableNamingConvention="structured",
        numberOfEventIndicators="0"
    )
    
    push!(root, h.ModelExchange(
        modelIdentifier=name,
        canGetAndSetFMUstate="true",
        canSerializeFMUstate="false",
        providesDirectionalDerivative="false"
    ))
    
    vars = h.ModelVariables()
    
    output_indices = Int[]
    derivative_indices = Int[]
    
    current_index = 1
    
    # 1. States x
    for i in 1:nx
        push!(vars, fmi2_real_variable("x[$i]", i, "local"; initial="exact", start=0.0))
        current_index += 1
    end
    
    # 2. Derivatives der(x)
    for i in 1:nx
        push!(vars, fmi2_real_variable("der(x[$i])", nx + i, "local"; derivative=i, initial="calculated"))
        push!(derivative_indices, nx + i)
        current_index += 1
    end
    
    # 3. Rigid body inputs
    vr_offset = 2 * nx
    for i in 1:3
        push!(vars, fmi2_real_variable("vb[$i]", vr_offset + i, "input"; start=0.0))
        current_index += 1
    end
    vr_offset += 3
    for i in 1:3
        push!(vars, fmi2_real_variable("wb[$i]", vr_offset + i, "input"; start=0.0))
        current_index += 1
    end
    vr_offset += 3
    for i in 1:3
        push!(vars, fmi2_real_variable("ab[$i]", vr_offset + i, "input"; start=0.0))
        current_index += 1
    end
    vr_offset += 3
    for i in 1:3
        push!(vars, fmi2_real_variable("dwb[$i]", vr_offset + i, "input"; start=0.0))
        current_index += 1
    end
    vr_offset += 3
    push!(vars, fmi2_real_variable("rho", vr_offset + 1, "input"; start=1.225))
    current_index += 1
    
    # 4. Control surfaces
    vr_offset += 1
    for k in 1:nCtrl
        cs_name = csNames[k]
        push!(vars, fmi2_real_variable(cs_name * "_delta", vr_offset + 1, "input"; start=0.0))
        push!(vars, fmi2_real_variable(cs_name * "_delta_dot", vr_offset + 2, "input"; start=0.0))
        push!(vars, fmi2_real_variable(cs_name * "_delta_ddot", vr_offset + 3, "input"; start=0.0))
        vr_offset += 3
        current_index += 3
    end
    
    # 5. Structural nodes
    if nVert > 0
        for i in 1:(3*nVert)
            push!(vars, fmi2_real_variable("rs[$i]", vr_offset + i, "input"; start=0.0))
            current_index += 1
        end
        vr_offset += 3 * nVert
        for i in 1:(3*nVert)
            push!(vars, fmi2_real_variable("vs[$i]", vr_offset + i, "input"; start=0.0))
            current_index += 1
        end
        vr_offset += 3 * nVert
        for i in 1:(3*nVert)
            push!(vars, fmi2_real_variable("as[$i]", vr_offset + i, "input"; start=0.0))
            current_index += 1
        end
        vr_offset += 3 * nVert
    end
    
    # 6. Outputs
    for i in 1:3
        push!(vars, fmi2_real_variable("Fb[$i]", vr_offset + i, "output"; initial="calculated"))
        push!(output_indices, current_index)
        current_index += 1
    end
    vr_offset += 3
    for i in 1:3
        push!(vars, fmi2_real_variable("Mb[$i]", vr_offset + i, "output"; initial="calculated"))
        push!(output_indices, current_index)
        current_index += 1
    end
    vr_offset += 3
    for k in 1:nmp
        mp_name = mpNames[k]
        for i in 1:6
            push!(vars, fmi2_real_variable(mp_name * "_load[$i]", vr_offset + i, "output"; initial="calculated"))
            push!(output_indices, current_index)
            current_index += 1
        end
        vr_offset += 6
    end
    
    push!(root, vars)
    
    # ModelStructure
    struct_node = h.ModelStructure()
    
    outputs_node = h.Outputs()
    for idx in output_indices
        push!(outputs_node, h.Unknown(index=string(idx)))
    end
    push!(struct_node, outputs_node)
    
    derivatives_node = h.Derivatives()
    for idx in derivative_indices
        push!(derivatives_node, h.Unknown(index=string(idx)))
    end
    push!(struct_node, derivatives_node)
    
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
