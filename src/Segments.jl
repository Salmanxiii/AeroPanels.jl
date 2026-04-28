
############################## Mesh Segments  #########################################

struct SegmentProperties{T}
    mid::Vector{Point3{T}}
    r::Vector{Point3{T}}
    aic3Ring::Matrix{Point3{T}}
    aic3Wake::Matrix{Point3{T}}
    panelToSeg::SparseMatrixCSC{Int, Int} 
end

function SegmentProperties{T}(nSegments::Int, nPanels::Int, nWakePanels::Int, panelToSeg::SparseMatrixCSC{Int, Int}) where T
    SegmentProperties(
    Vector{Point3{T}}(undef, nSegments),
    Vector{Point3{T}}(undef, nSegments),
    Matrix{Point3{T}}(undef, nSegments, nPanels),
    Matrix{Point3{T}}(undef, nSegments, nWakePanels),
    panelToSeg)
end

function ProcessSegments(meshRing::GeometryBasics.Mesh{3, T}, sizesRing, meshWake::GeometryBasics.Mesh{3, T}, sizesWake, symmXZ) where T
    nss, ncs = TotalSegments(sizesRing)
    nr = sizesRing.totalPanels
    nw = sizesWake.totalPanels
    segmentSpan = SegmentProperties{T}(nss, nr, nw, BuildSpanwiseCirculationOperator(sizesRing))
    segmentChord = SegmentProperties{T}(ncs, nr, nw, BuildChordwiseCirculationOperator(sizesRing))
    ProcessSegments!(segmentSpan, segmentChord, meshRing, sizesRing, meshWake, sizesWake, symmXZ)
    return segmentSpan, segmentChord
end

function ProcessSegments!(segmentSpan, segmentChord, meshRing, sizesRing, meshWake, sizesWake, symmXZ)
#function ProcessSegments!(segmentSpan, segmentChord, vertices, sizesRing)
    vertices = coordinates(meshRing)
    # loop over Spanwise Segments
    for (s, nc, ns) in sizesRing
        @batch for i in 1:nc+1
            for j in 1:ns
                r1 = vertices[VertexIndex(s, i, j, sizesRing)]
                r2 = vertices[VertexIndex(s, i, j+1, sizesRing)]
                m = SpanSegmentIndex(s, i, j, sizesRing)
                segmentSpan.mid[m] = (r1 + r2) / 2
                segmentSpan.r[m] = r2 - r1
            end
        end
    end

    # loop over Chordwise Segments
    for (s, nc, ns) in sizesRing
        @batch for i in 1:nc
            for j in 1:ns+1
                r1 = vertices[VertexIndex(s, i, j, sizesRing)]
                r2 = vertices[VertexIndex(s, i+1, j, sizesRing)]
                m = ChordSegmentIndex(s, i, j, sizesRing)
                segmentChord.mid[m] = (r1 + r2) / 2
                segmentChord.r[m] = r2 - r1
            end    
        end
    end

    Influence3!(segmentSpan.aic3Ring, segmentSpan.mid, meshRing, symmXZ)
    Influence3!(segmentSpan.aic3Wake, segmentSpan.mid, meshWake, symmXZ)
    Influence3!(segmentChord.aic3Ring, segmentChord.mid, meshRing, symmXZ)
    Influence3!(segmentChord.aic3Wake, segmentChord.mid, meshWake, symmXZ)
end

function SegmentCirculation(Γp::Vector{T}, sizes::Sizes) where T
    Γs = Vector{T}(undef, sizes.totalSpanSegments)
    Γc = Vector{T}(undef, sizes.totalChordSegments)
    SegmentCirculation!(Γs, Γc, Γp, sizes)
    return Γs, Γc

end

function SegmentCirculation!(Γs, Γc, Γp, sizes::Sizes)
    SpanSegmentCirculation!(Γs, Γp, sizes)
    ChordSegmentCirculation!(Γc, Γp, sizes)
end

function SpanSegmentCirculation!(Γs, Γp, sizes::Sizes)
    for (s, nc, ns) in sizes
        @batch for i in 1:nc+1
            for j in 1:ns
                m = SpanSegmentIndex(s, i, j, sizes)
                Γs[m] = i == 1 ? Γp[PanelIndex(s, i, j, sizes)] :
                        i == nc+1 ? 0 :
                        Γp[PanelIndex(s, i, j, sizes)] - Γp[PanelIndex(s, i-1, j, sizes)]
            end
        end
    end
end

function ChordSegmentCirculation!(Γc, Γp, sizes::Sizes)
    for (s, nc, ns) in sizes
        @batch for i in 1:nc
            for j in 1:ns+1
                m = ChordSegmentIndex(s, i, j, sizes)
                Γc[m] = j == 1 ? Γp[PanelIndex(s, i, j, sizes)] :
                        j == ns+1 ? -Γp[PanelIndex(s, i, ns, sizes)] :
                        Γp[PanelIndex(s, i, j, sizes)] - Γp[PanelIndex(s, i, j-1, sizes)]
            end
        end
    end
end

"""
Returns a sparse operator that maps panel circulations to spanwise segment circulations
Γs = operator * Γp
"""
function BuildSpanwiseCirculationOperator(sizes::Sizes)
    n_pos = sizes.totalPanels
    n_neg = sum((nc - 1) * ns for (_, nc, ns) in sizes)
    pos_from, pos_to = Vector{Int}(undef, n_pos), Vector{Int}(undef, n_pos)
    neg_from, neg_to = Vector{Int}(undef, n_neg), Vector{Int}(undef, n_neg)
    indp, indn = 1, 1
    for (s, nc, ns) in sizes, j in 1:ns
        for i in 1:nc
            p_idx = PanelIndex(s, i, j, sizes)
            # Term 1: + Γp[i] -> Maps Panel i to Segment i
            pos_from[indp], pos_to[indp] = p_idx, SpanSegmentIndex(s, i, j, sizes)
            indp += 1
            # Term 2: - Γp[i] -> Maps Panel i to Segment i+1 (effectively -Γp[i-1] at next row)
            if i < nc
                neg_from[indn], neg_to[indn] = p_idx, SpanSegmentIndex(s, i+1, j, sizes)
                indn += 1
            end
        end
    end
    M_pos = SelectionOperator(pos_from, sizes.totalPanels, pos_to, sizes.totalSpanSegments, Int)
    M_neg = SelectionOperator(neg_from, sizes.totalPanels, neg_to, sizes.totalSpanSegments, Int)
    return M_pos - M_neg
end

function BuildChordwiseCirculationOperator(sizes::Sizes)
    n_entries = sizes.totalPanels
    pos_from, pos_to = Vector{Int}(undef, n_entries), Vector{Int}(undef, n_entries)
    neg_from, neg_to = Vector{Int}(undef, n_entries), Vector{Int}(undef, n_entries)

    indp, indn = 1, 1
    for (s, nc, ns) in sizes, i in 1:nc, j in 1:ns
        p_idx = PanelIndex(s, i, j, sizes)
        # Term 1: + Γp[i, j] -> Maps Panel j to Segment j
        pos_from[indp], pos_to[indp] = p_idx, ChordSegmentIndex(s, i, j, sizes)
        indp += 1
        # Term 2: - Γp[i, j] -> Maps Panel j to Segment j+1 (effectively -Γp[i, j-1] at next row)
        neg_from[indn], neg_to[indn] = p_idx, ChordSegmentIndex(s, i, j+1, sizes)
        indn += 1
    end
    M_pos = SelectionOperator(pos_from, sizes.totalPanels, pos_to, sizes.totalChordSegments, Int)
    M_neg = SelectionOperator(neg_from, sizes.totalPanels, neg_to, sizes.totalChordSegments, Int)

    return M_pos - M_neg
end


"""
Calculate the velocity at each segment due to the ring and wake influences
"""
function SegmentInducedVelocity(Γr, Γw, sp::SegmentProperties)
    # Non-allocating version of equation v = AIC3r * Γr + AIC3w * Γw
    V = sp.aic3Ring * Γr
    mul!(V, sp.aic3Wake, Γw, 1.0, 1.0)
end
"""
Calculate the force on each segment due to the circulation and airflow velocity at the segment
"""
function SegmentForce(Γr, Γw, vb, ρ, sp::SegmentProperties)
    Γs = sp.panelToSeg * Γr
    Fs = SegmentInducedVelocity(Γr, Γw, sp) # calculate velocity induced at each segment
    # Non-allocating version of Fs = ρ * Γs .* cross.(Vs, rs)
    @batch for i in 1:length(Fs)
        Fs[i] = ρ * Γs[i] * cross(Fs[i]+vb, sp.r[i]) 
    end
    return Fs
end
function SegmentForce(Γr, vb, ρ, sp::SegmentProperties)
    Γs = sp.panelToSeg * Γr
    Fs = [ρ * Γs[i] * cross(vb, sp.r[i]) for i in 1:length(Γs)]
    return Fs
end
