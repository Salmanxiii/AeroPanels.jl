############################## Indexing Struct ####################################

"""
$(TYPEDEF)

Struct to store the sizes and pre-calculated indices of the mesh elements.

$(TYPEDFIELDS)
"""
struct Sizes
    sizes::Vector{Tuple{Int, Int}}
    nSurfaces::Int
    totalPanels::Int
    totalVertices::Int
    totalSpanSegments::Int
    totalChordSegments::Int
    panelIndices::Vector{Matrix{Int}}
    vertexIndices::Vector{Matrix{Int}}
    spanSegmentIndices::Vector{Matrix{Int}}
    chordSegmentIndices::Vector{Matrix{Int}}
end

function Base.iterate(s::Sizes, state=1)
    state > length(s) ? nothing : ((state, s.sizes[state][1], s.sizes[state][2]), state + 1)
end

function Base.length(s::Sizes)
    return s.nSurfaces
end

# return a vector of matrix of indices for each surface in the mesh
# used to index into faces and vertices of mesh
function IndicesMatrix(sizesVec::Vector{Tuple{Int, Int}})
    indices = [Matrix{Int}(undef, nc, ns) for (nc, ns) in sizesVec]
    adder = 0
    for (s,(nc,ns)) in enumerate(sizesVec)
        for i in 1:nc, j in 1:ns
            indices[s][i, j] = adder + ns * (i - 1) + j
        end
        adder += nc*ns
    end
    return indices
end

# Constructor for the Sizes struct
function Sizes(sizesVec::Vector{Tuple{Int, Int}})
    totalPanels = sum(nc * ns for (nc, ns) in sizesVec)
    totalVertices = sum((nc + 1) * (ns + 1) for (nc, ns) in sizesVec)
    totalSpanSegments = sum((nc+1)*ns for (nc, ns) in sizesVec)
    totalChordSegments = sum(nc*(ns+1) for (nc, ns) in sizesVec)
    panelIndices = IndicesMatrix(sizesVec)
    vertexIndices = IndicesMatrix([(nc+1, ns+1) for (nc, ns) in sizesVec])
    spanSegmentIndices = IndicesMatrix([(nc+1, ns) for (nc, ns) in sizesVec])
    chordSegmentIndices = IndicesMatrix([(nc, ns+1) for (nc, ns) in sizesVec])
    return Sizes(sizesVec, length(sizesVec), totalPanels, totalVertices, totalSpanSegments,
     totalChordSegments, panelIndices, vertexIndices, spanSegmentIndices, chordSegmentIndices)
end

"""
    PanelIndex(s::Int, i::Int, j::Int, sizes::Sizes)

Return the global index of the panel `(i, j)` on surface `s`.
"""
PanelIndex(s::Int, i::Int, j::Int, sizes::Sizes) = sizes.panelIndices[s][i, j]

"""
    TEPanelIndex(sizes::Sizes)

Return a vector of global indices for all panels along the trailing edges.
"""
TEPanelIndex(sizes::Sizes) = [PanelIndex(s, nc, j, sizes) for (s, nc, ns) in sizes for j in 1:ns]

"""
    LEPanelIndex(sizes::Sizes)

Return a vector of global indices for all panels along the leading edges.
"""
LEPanelIndex(sizes::Sizes) = [PanelIndex(s, 1, j, sizes)  for (s, nc, ns) in sizes for j in 1:ns]

NonKuttaPanelIndex(sizes::Sizes) = [PanelIndex(s, i, j, sizes)  for (s, nc, ns) in sizes for i in 2:nc for j in 1:ns]

VertexIndex(s::Int, i::Int, j::Int, sizes::Sizes) = sizes.vertexIndices[s][i, j]
TEVertexIndex(sizes::Sizes) = [VertexIndex(s, nc+1, j, sizes) for (s, nc, ns) in sizes for j in 1:ns+1]
LEVertexIndex(sizes::Sizes) = [VertexIndex(s, 1, j, sizes)    for (s, nc, ns) in sizes for j in 1:ns+1]

SpanSegmentIndex(s::Int, i::Int, j::Int, sizes::Sizes) = sizes.spanSegmentIndices[s][i, j]
ChordSegmentIndex(s::Int, i::Int, j::Int, sizes::Sizes) = sizes.chordSegmentIndices[s][i, j]
TotalSegments(sizes::Sizes) = sizes.totalSpanSegments, sizes.totalChordSegments

function SelectionOperator(fromIndices::AbstractVector{Int}, fromArraySize::Int, 
    toIndices::AbstractVector{Int}, toArraySize::Int, T::Type=Float64)
   n = length(toIndices)
   @assert length(toIndices) == length(fromIndices)
   op = sparse(toIndices, fromIndices, ones(T, n), toArraySize, fromArraySize) 
   return op
end