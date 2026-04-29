struct SegmentProperties{T}
    mid::Vector{Point3{T}}
    r::Vector{Point3{T}}
    aic3Ring::Matrix{Point3{T}}
    aic3Wake::Matrix{Point3{T}}
    indPosSeg::Vector{Int}
    indPosPanel::Vector{Int}
    indNegSeg::Vector{Int}
    indNegPanel::Vector{Int}
    nSpanSegments::Int
    nChordSegments::Int
    nTotalSegments::Int
end

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

function BuildSegmentIndices(sizes::Sizes)
    nss, ncs = TotalSegments(sizes)
    pos_seg1, pos_pan1, neg_seg1, neg_pan1 = BuildSpanwiseIndices(sizes)
    pos_seg2, pos_pan2, neg_seg2, neg_pan2 = BuildChordwiseIndices(sizes)
    pos_seg2 .+= nss
    neg_seg2 .+= nss
    return [pos_seg1; pos_seg2], [pos_pan1; pos_pan2], [neg_seg1; neg_seg2], [neg_pan1; neg_pan2]
end

function ProcessSegments(meshRing::GeometryBasics.Mesh{3, T}, sizesRing, meshWake::GeometryBasics.Mesh{3, T}, sizesWake, symmXZ) where T
    nss, ncs = TotalSegments(sizesRing)
    nr = sizesRing.totalPanels
    nw = sizesWake.totalPanels
    nts = nss + ncs # number of total segments in body

    pos_seg, pos_pan, neg_seg, neg_pan = BuildSegmentIndices(sizesRing)

    segmentProps = SegmentProperties(
    Vector{Point3{T}}(undef, nts),
    Vector{Point3{T}}(undef, nts),
    Matrix{Point3{T}}(undef, nts, nr),
    Matrix{Point3{T}}(undef, nts, nw),
    pos_seg, pos_pan, neg_seg, neg_pan, nss, ncs, nts)

    UpdateSegmentProperties!(segmentProps, meshRing, sizesRing)
    UpdateSegmentAIC!(segmentProps, meshRing, meshWake, symmXZ)
    return segmentProps
end

function UpdateSegmentProperties!(segmentProps, meshRing, sizesRing)
    vertices = coordinates(meshRing)
    nss = segmentProps.nSpanSegments

    # loop over Spanwise Segments
    for (s, nc, ns) in sizesRing
        @batch for i in 1:nc+1
            for j in 1:ns
                r1 = vertices[VertexIndex(s, i, j, sizesRing)]
                r2 = vertices[VertexIndex(s, i, j+1, sizesRing)]
                m = SpanSegmentIndex(s, i, j, sizesRing)
                segmentProps.mid[m] = (r1 + r2) / 2
                segmentProps.r[m] = r2 - r1
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
                segmentProps.mid[m+nss] = (r1 + r2) / 2
                segmentProps.r[m+nss] = r2 - r1
            end    
        end
    end
    return nothing
end

function UpdateSegmentAIC!(segmentProps, meshRing, meshWake, symmXZ)
    Influence3!(segmentProps.aic3Ring, segmentProps.mid, meshRing, symmXZ)
    Influence3!(segmentProps.aic3Wake, segmentProps.mid, meshWake, symmXZ)
    return nothing
end

function SegmentCirculation!(Γseg, Γp, sp::SegmentProperties)
    Γseg .= 0
    Γseg[sp.indPosSeg] .= @view Γp[sp.indPosPanel]
    Γseg[sp.indNegSeg] .-= @view Γp[sp.indNegPanel]
    return Γseg
end
SegmentCirculation(Γp::AbstractArray{T}, sp::SegmentProperties) where T = SegmentCirculation!(zeros(T, size(sp.aic3Ring,1)), Γp, sp)
