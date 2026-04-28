
# 1. New struct definition (Iteration 2 as requested)
struct SegmentProperties{T}
    mid::Vector{Point3{T}}
    r::Vector{Point3{T}}
    aic3Ring::Matrix{Point3{T}}
    aic3Wake::Matrix{Point3{T}}
    indPosSeg::Vector{Int}
    indPosPanel::Vector{Int}
    indNegSeg::Vector{Int}
    indNegPanel::Vector{Int}
end

# 2. Population logic for indices
function BuildSpanwiseIndices(sizes::Sizes)
    n_pos = sizes.totalPanels
    n_neg = sum((nc - 1) * ns for (_, nc, ns) in sizes)
    pos_pan, pos_seg = Vector{Int}(undef, n_pos), Vector{Int}(undef, n_pos)
    neg_pan, neg_seg = Vector{Int}(undef, n_neg), Vector{Int}(undef, n_neg)
    indp, indn = 1, 1
    for (s, nc, ns) in sizes, j in 1:ns
        for i in 1:nc
            p_idx = PanelIndex(s, i, j, sizes)
            # Term 1: + Γp[i] -> Maps Panel i to Segment i
            pos_pan[indp], pos_seg[indp] = p_idx, SpanSegmentIndex(s, i, j, sizes)
            indp += 1
            # Term 2: - Γp[i] -> Maps Panel i to Segment i+1 (effectively -Γp[i-1] at next row)
            if i < nc
                neg_pan[indn], neg_seg[indn] = p_idx, SpanSegmentIndex(s, i+1, j, sizes)
                indn += 1
            end
        end
    end
    return pos_seg, pos_pan, neg_seg, neg_pan
end

function BuildChordwiseIndices(sizes::Sizes)
    n_entries = sizes.totalPanels
    pos_pan, pos_seg = Vector{Int}(undef, n_entries), Vector{Int}(undef, n_entries)
    neg_pan, neg_seg = Vector{Int}(undef, n_entries), Vector{Int}(undef, n_entries)

    indp, indn = 1, 1
    for (s, nc, ns) in sizes, i in 1:nc, j in 1:ns
        p_idx = PanelIndex(s, i, j, sizes)
        # Term 1: + Γp[i, j] -> Maps Panel j to Segment j
        pos_pan[indp], pos_seg[indp] = p_idx, ChordSegmentIndex(s, i, j, sizes)
        indp += 1
        # Term 2: - Γp[i, j] -> Maps Panel j to Segment j+1 (effectively -Γp[i, j-1] at next row)
        neg_pan[indn], neg_seg[indn] = p_idx, ChordSegmentIndex(s, i, j+1, sizes)
        indn += 1
    end
    return pos_seg, pos_pan, neg_seg, neg_pan
end

function ProcessSegments(meshRing::GeometryBasics.Mesh{3, T}, sizesRing, meshWake::GeometryBasics.Mesh{3, T}, sizesWake, symmXZ) where T
    nss, ncs = TotalSegments(sizesRing)
    nr = sizesRing.totalPanels
    nw = sizesWake.totalPanels

    pos_seg, pos_pan, neg_seg, neg_pan = BuildSpanwiseIndices(sizesRing)
    segmentSpan = SegmentProperties(
    Vector{Point3{T}}(undef, nss),
    Vector{Point3{T}}(undef, nss),
    Matrix{Point3{T}}(undef, nss, nr),
    Matrix{Point3{T}}(undef, nss, nw),
    pos_seg, pos_pan, neg_seg, neg_pan)

    pos_seg, pos_pan, neg_seg, neg_pan = BuildChordwiseIndices(sizesRing)
    segmentChord = SegmentProperties(
    Vector{Point3{T}}(undef, ncs),
    Vector{Point3{T}}(undef, ncs),
    Matrix{Point3{T}}(undef, ncs, nr),
    Matrix{Point3{T}}(undef, ncs, nw),
    pos_seg, pos_pan, neg_seg, neg_pan)

    UpdateSegmentProperties!(segmentSpan, segmentChord, meshRing, sizesRing)
    UpdateSegmentAIC!(segmentSpan, segmentChord, meshRing, meshWake, symmXZ)
    return segmentSpan, segmentChord
end

function UpdateSegmentProperties!(segmentSpan, segmentChord, meshRing, sizesRing)
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
    return nothing
end

function UpdateSegmentAIC!(segmentSpan, segmentChord, meshRing, meshWake, symmXZ)
    Influence3!(segmentSpan.aic3Ring, segmentSpan.mid, meshRing, symmXZ)
    Influence3!(segmentSpan.aic3Wake, segmentSpan.mid, meshWake, symmXZ)
    Influence3!(segmentChord.aic3Ring, segmentChord.mid, meshRing, symmXZ)
    Influence3!(segmentChord.aic3Wake, segmentChord.mid, meshWake, symmXZ)
    return nothing
end

# 3. Calculation logic using the new indexing scheme
function SegmentCirculation!(Γseg, Γp, sp::SegmentProperties)
    Γseg .= 0
    Γseg[sp.indPosSeg] .= @view Γp[sp.indPosPanel]
    Γseg[sp.indNegSeg] .-= @view Γp[sp.indNegPanel]
    return Γseg
end
SegmentCirculation(Γp::AbstractArray{T}, sp::SegmentProperties) where T = SegmentCirculation!(zeros(T, size(sp.aic3Ring,1)), Γp, sp)



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
    Γs = SegmentCirculation(Γr, sp)
    Fs = SegmentInducedVelocity(Γr, Γw, sp) # calculate velocity induced at each segment
    # Non-allocating version of Fs = ρ * Γs .* cross.(Vs, rs)
    @batch for i in 1:length(Fs)
        Fs[i] = ρ * Γs[i] * cross(Fs[i]+vb, sp.r[i]) 
    end
    return Fs
end
function SegmentForce(Γr, vb, ρ, sp::SegmentProperties)
    Γs = SegmentCirculation!(Γs, Γr, sp)
    Fs = [ρ * Γs[i] * cross(vb, sp.r[i]) for i in 1:length(Γs)]
    return Fs
end
