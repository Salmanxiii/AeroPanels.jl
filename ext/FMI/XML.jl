module FMI_XML

using XML

export generate_xml_node

function generate_xml_node(name, nx, version)
    h = XML.h
    root = h.fmiModelDescription(
        fmiVersion="3.0", modelName=name,
        instantiationToken="{$(bytes2hex(rand(UInt8, 16)))}",
        generationTool="AeroPanels.jl", version=version
    )
    push!(root, h.ModelExchange(modelIdentifier=name))
    push!(root, h.CoSimulation(modelIdentifier=name))
    vars = h.ModelVariables()
    push!(vars, h.Float64(name="time", valueReference="0", causality="independent", variability="continuous"))
    
    x_var = h.Float64(name="x", valueReference="1", causality="local", variability="continuous", initial="exact")
    push!(x_var, h.Dimension(start=string(nx)))
    x_var["start"] = join(fill(0.0, nx), " ")
    push!(vars, x_var)

    u_var = h.Float64(name="u", valueReference="2", causality="input", variability="continuous")
    push!(u_var, h.Dimension(start="13"))
    u_var["start"] = join([fill(0.0, 12)..., 1.225], " ")
    push!(vars, u_var)

    y_var = h.Float64(name="y", valueReference="3", causality="output", variability="continuous", initial="calculated")
    push!(y_var, h.Dimension(start="6"))
    push!(vars, y_var)

    der_x = h.Float64(name="der(x)", valueReference="4", causality="local", variability="continuous", derivative="1", initial="calculated")
    push!(der_x, h.Dimension(start=string(nx)))
    push!(vars, der_x)

    push!(root, vars)
    struct_node = h.ModelStructure()
    push!(struct_node, h.Output(valueReference="3"))
    push!(struct_node, h.ContinuousStateDerivative(valueReference="4"))
    push!(root, struct_node)
    return root
end

end
