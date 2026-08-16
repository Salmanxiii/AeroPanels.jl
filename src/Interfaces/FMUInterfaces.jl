abstract type AbstractFMUModel{T} <: AeroModel end
abstract type AbstractFMUArray{T} <: AbstractVector{T} end

# Mandatory interface contract for all FMU-compatible models
GetTotalSize(model::AbstractFMUModel) = error("GetTotalSize not implemented for $(typeof(model))")
CreateFMUArray(model::AbstractFMUModel) = error("CreateFMUArray not implemented for $(typeof(model))")
GetFMULayout(model::AbstractFMUModel) = error("GetFMULayout not implemented for $(typeof(model))")
NumberOfStates(model::AbstractFMUModel) = error("NumberOfStates not implemented for $(typeof(model))")

AllocateFMUCaches(model::AbstractFMUModel, ::Type{T}) where T = error("AllocateFMUCaches not implemented for $(typeof(model))")
AllocateFMUCaches(model::AbstractFMUModel) = AllocateFMUCaches(model, Float64)

InitializeFMU!(array::AbstractFMUArray, cache, model::AbstractFMUModel, t::Float64; start_from_trim::Bool=false) = error("InitializeFMU! not implemented for $(typeof(model))")
InitializeSteadyState!(array::AbstractFMUArray, cache, model::AbstractFMUModel, t::Float64) = error("InitializeSteadyState! not implemented for $(typeof(model))")
EvaluateDerivatives!(array::AbstractFMUArray, cache, model::AbstractFMUModel, t::Float64) = error("EvaluateDerivatives! not implemented for $(typeof(model))")
function EvaluateDerivatives!(du::AbstractVector, array::AbstractFMUArray, cache, model::AbstractFMUModel, t::Float64)
    EvaluateDerivatives!(array, cache, model, t)
    nStates = NumberOfStates(model)
    if nStates > 0
        du .= view(array.data, (nStates + 1):(2 * nStates))
    end
    return nothing
end

EvaluateOutputs!(array::AbstractFMUArray, cache, model::AbstractFMUModel, t::Float64) = error("EvaluateOutputs! not implemented for $(typeof(model))")


ManualDirectionalDerivative(FMUmodel::AbstractFMUModel, FMUarray::AbstractFMUArray, outputRefs::AbstractVector{<:Integer}, inputRefs::AbstractVector{<:Integer}, inputSeed::AbstractVector; cache=nothing) = nothing

function GetDirectionalDerivative!(
    d_dx::AbstractVector{T},
    outputRefs::AbstractVector{<:Integer},
    inputRefs::AbstractVector{<:Integer},
    inputSeed::AbstractVector,
    FMUarray::AbstractFMUArray{ForwardDiff.Dual{Tag, T, 1}},
    FMUmodel::AbstractFMUModel,
    cache,
    t::Float64=0.0
) where {Tag, T}
    manual_deriv = ManualDirectionalDerivative(FMUmodel, FMUarray, outputRefs, inputRefs, inputSeed; cache=cache)
    if !isnothing(manual_deriv)
        d_dx .= manual_deriv
        return nothing
    end

    layout = GetFMULayout(FMUmodel)
    nStates = NumberOfStates(FMUmodel)
    out_start = layout.OutputsStartIndex
    # Flush all array partials before seeding
    @inbounds for idx in 1:length(FMUarray.data)
        FMUarray.data[idx] = ForwardDiff.Dual{Tag, T, 1}(ForwardDiff.value(FMUarray.data[idx]), ForwardDiff.Partials((zero(T),)))
    end

    # Seed input/state partials
    @inbounds for i in 1:length(inputRefs)
        idx = inputRefs[i]
        FMUarray.data[idx] = ForwardDiff.Dual{Tag, T, 1}(ForwardDiff.value(FMUarray.data[idx]), ForwardDiff.Partials((T(inputSeed[i]),)))
    end

    # Determine required evaluations
    has_derivatives = false
    has_outputs = false
    out_start = layout.OutputsStartIndex
    @inbounds for ref in outputRefs
        if ref < out_start
            has_derivatives = true
        else
            has_outputs = true
        end
    end

    if has_derivatives
        EvaluateDerivatives!(FMUarray, cache, FMUmodel, t)
    end

    if has_outputs
        EvaluateOutputs!(FMUarray, cache, FMUmodel, t)
    end

    # Extract output/derivative partials
    @inbounds for k in 1:length(outputRefs)
        ref = outputRefs[k]
        d_dx[k] = ForwardDiff.partials(FMUarray.data[ref], 1)
    end

    return nothing
end

function GetFMISparsity(model::AbstractFMUModel)
    error("GetFMISparsity is not yet implemented for $(typeof(model)). Automated sparsity tracing via SparseConnectivityTracer is a work-in-progress ToDo.")
end

Base.size(a::AbstractFMUArray) = size(a.data)
Base.getindex(a::AbstractFMUArray, i::Int) = a.data[i]
Base.setindex!(a::AbstractFMUArray, v, i::Int) = (a.data[i] = v)
Base.pointer(a::AbstractFMUArray) = pointer(a.data)
Base.IndexStyle(::Type{<:AbstractFMUArray}) = IndexLinear()

struct FMULayout
    ArrayType::DataType
    States::Tuple
    Derivatives::Tuple
    Inputs::Tuple
    Outputs::Tuple
    OutputsStartIndex::Int
    OutputStates::Tuple
    OutputDerivatives::Tuple
end

# Backwards compatible constructor
function FMULayout(ArrayType::DataType, States::Tuple, Derivatives::Tuple, Inputs::Tuple, Outputs::Tuple, OutputsStartIndex::Int)
    return FMULayout(ArrayType, States, Derivatives, Inputs, Outputs, OutputsStartIndex, (), ())
end

macro fmu_model(expr)
    if !(expr.head == :struct)
        error("Expected a struct definition")
    end
    
    # 1. Parse name and type parameters
    sig = expr.args[2]
    struct_head = (sig isa Expr && sig.head == :<:) ? sig.args[1] : sig
    
    local name_sym
    local type_params = []
    
    if struct_head isa Symbol
        name_sym = struct_head
    elseif struct_head isa Expr && struct_head.head == :curly
        name_sym = struct_head.args[1]
        type_params = struct_head.args[2:end]
    end
    
    local type_params_clean = [p isa Expr && p.head == :<: ? p.args[1] : p for p in type_params]
    local T_sym_with_bounds = isempty(type_params) ? :T : type_params[end]
    local T_sym = (T_sym_with_bounds isa Expr && T_sym_with_bounds.head == :<:) ? T_sym_with_bounds.args[1] : T_sym_with_bounds
    
    # 2. Parse fields and separate layout
    body = expr.args[3]
    physics_fields = []
    
    states_def = []       # Tuple of (fname, fval, is_output)
    derivatives_def = []  # Tuple of (fname, fval, is_output)
    inputs_def = []
    outputs_def = []
    
    current_section = :none
    
    for arg in body.args
        if arg isa Expr && arg.head == :vect && length(arg.args) == 1 && arg.args[1] isa Symbol
            sym = arg.args[1]
            if sym == :States
                current_section = :states
            elseif sym == :Derivatives
                current_section = :derivatives
            elseif sym == :Inputs
                current_section = :inputs
            elseif sym == :Outputs
                current_section = :outputs
            else
                push!(physics_fields, arg)
            end
        elseif arg isa Expr && (arg.head == :(=) || arg.head == :kw)
            fname = arg.args[1]
            fval_raw = arg.args[2]
            
            # Helper to parse (size=..., output=true) or plain size
            fval = fval_raw
            is_output = false
            if fval_raw isa Expr && fval_raw.head == :tuple
                for kw in fval_raw.args
                    if kw isa Expr && (kw.head == :(=) || kw.head == :kw)
                        key, val = kw.args[1], kw.args[2]
                        if key == :size
                            fval = val
                        elseif key == :output
                            is_output = (val == true)
                        end
                    end
                end
            end
            
            if current_section == :states
                push!(states_def, (fname, fval, is_output))
            elseif current_section == :derivatives
                push!(derivatives_def, (fname, fval, is_output))
            elseif current_section == :inputs
                push!(inputs_def, (fname, fval, false))
            elseif current_section == :outputs
                push!(outputs_def, (fname, fval, true))
            else
                push!(physics_fields, arg)
            end
        else
            if current_section == :none
                push!(physics_fields, arg)
            end
        end
    end
    
    # Replace body args with just physics fields
    new_body = Expr(:block, physics_fields...)
    new_struct = Expr(:struct, expr.args[1], expr.args[2], new_body)
    
    # Generate companion array name
    array_name = Symbol(string(name_sym) * "FMUArray")
    
    all_defs = vcat(states_def, derivatives_def, inputs_def, outputs_def)
    
    # Generate array struct
    array_fields = Expr(:block, :(data::Vector{V}))
    for (fname, _, _) in all_defs
        push!(array_fields.args, :($fname::SubArray{V, 1, Vector{V}, Tuple{UnitRange{Int}}, true}))
    end
    
    # Generate CreateFMUArray
    total_len_expr = :(0)
    for (_, fval, _) in all_defs
        if fval isa Int
            total_len_expr = :($total_len_expr + $fval)
        else
            total_len_expr = :($total_len_expr + ($fval)(m))
        end
    end
    
    ctor_body_with_data = Expr(:block)
    push!(ctor_body_with_data.args, :(current_idx = 1))
    
    for (fname, fval, _) in all_defs
        push!(ctor_body_with_data.args, quote
            sz = $(fval) isa Int ? $(fval) : ($(fval))(m)
            end_idx = current_idx + sz - 1
            $fname = view(data, current_idx:end_idx)
            current_idx = end_idx + 1
        end)
    end
    
    return_args = Any[:data]
    for (fname, _, _) in all_defs
        push!(return_args, fname)
    end
    
    array_type_expr = :($array_name{V})
    array_decl_expr = :($array_name{V<:Real})
    model_type_expr = isempty(type_params_clean) ? name_sym : :($name_sym{$(type_params_clean...)})
    
    push!(ctor_body_with_data.args, :(return $array_type_expr($(return_args...))))
    
    # Layout metadata
    states_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,_) in states_def]...)
    derivatives_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,_) in derivatives_def]...)
    inputs_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,_) in inputs_def]...)
    outputs_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,_) in outputs_def]...)
    
    out_states_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,out) in states_def if out]...)
    out_derivs_tup = Expr(:tuple, [QuoteNode(fname) for (fname,_,out) in derivatives_def if out]...)
    
    layout_array_type = :($array_name{Float64})
    layout_body = Expr(:block)
    push!(layout_body.args, :(out_start = 1))
    for (_, fval, _) in vcat(states_def, derivatives_def, inputs_def)
        push!(layout_body.args, :(out_start += $(fval) isa Int ? $(fval) : ($(fval))(m)))
    end
    push!(layout_body.args, :(return FMULayout($layout_array_type, $states_tup, $derivatives_tup, $inputs_tup, $outputs_tup, out_start, $out_states_tup, $out_derivs_tup)))

    # 4. Final generation with a single escape
    return esc(quote
        Base.@__doc__ $new_struct
        
        struct $array_decl_expr <: AbstractFMUArray{V}
            $array_fields
        end
        
        function AeroPanels.CreateFMUArray(m::$model_type_expr, data::AbstractVector{V}) where {$(type_params_clean...), V}
            $ctor_body_with_data
        end

        function AeroPanels.CreateFMUArray(::Type{V}, m::$model_type_expr) where {$(type_params_clean...), V}
            totalLen = $total_len_expr
            data = zeros(V, totalLen)
            return AeroPanels.CreateFMUArray(m, data)
        end

        function AeroPanels.CreateFMUArray(m::$model_type_expr) where {$(type_params_clean...)}
            return AeroPanels.CreateFMUArray($T_sym, m)
        end
        
        AeroPanels.GetTotalSize(m::$model_type_expr) where {$(type_params_clean...)} = length(AeroPanels.CreateFMUArray(m).data)
        AeroPanels.NumberOfStates(m::$model_type_expr) where {$(type_params_clean...)} = isempty(Int[$( [:( $(fval) isa Int ? $(fval) : ($(fval))(m) ) for (_, fval, _) in states_def]... )]) ? 0 : sum(Int[$( [:( $(fval) isa Int ? $(fval) : ($(fval))(m) ) for (_, fval, _) in states_def]... )])
        
        function AeroPanels.GetFMULayout(m::$model_type_expr) where {$(type_params_clean...)}
            $layout_body
        end
    end)
end
