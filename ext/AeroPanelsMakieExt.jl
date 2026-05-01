module AeroPanelsMakieExt

using GLMakie
using AeroPanels
using GeometryBasics
using StaticArrays

function disconnected_mesh(mesh, face_colors)
    fs = faces(mesh)
    vs = coordinates(mesh)
    new_vs = similar(vs, 0)
    new_fs = similar(fs, 0)
    new_cs = similar(face_colors, 0)
    
    idx = 1
    for (i, f) in enumerate(fs)
        for v in f
            push!(new_vs, vs[v])
            push!(new_cs, face_colors[i])
        end
        # Assuming QuadFace
        push!(new_fs, typeof(f)(idx, idx+1, idx+2, idx+3))
        idx += 4
    end
    return GeometryBasics.Mesh(new_vs, new_fs), new_cs
end

function AeroPanels.PlotModel(model::AeroModel; 
    plotWake=false, 
    Γp=nothing, 
    Γw=nothing, 
    plotForces=false, 
    sol=nothing, 
    forceScale=1.0)
    
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    # Plot panel circulation if provided, else wireframe
    if !isnothing(Γp)
        flat_mesh_p, v_colors_p = disconnected_mesh(model.mesh, Γp)
        mesh!(ax, flat_mesh_p, color=v_colors_p, colormap=:viridis)
        wireframe!(ax, model.mesh, color=(:black, 0.5), linewidth=0.5)
    else
        wireframe!(ax, model.mesh, color=:blue, linewidth=1.5)
    end

    if plotWake
        if !isnothing(Γw)
            flat_mesh_w, v_colors_w = disconnected_mesh(model.wakeMesh, Γw)
            mesh!(ax, flat_mesh_w, color=v_colors_w, colormap=:viridis)
            wireframe!(ax, model.wakeMesh, color=(:black, 0.5), linewidth=0.5)
        else
            wireframe!(ax, model.wakeMesh, color=:gray, linewidth=0.5)
        end
    end

    if plotForces && !isnothing(sol)
        locations = AeroPanels.AerodynamicLoadLocation(model)
        forces = sol.forceVecSeg
        dirs = [Point3f(f[1], f[2], f[3]) for f in forces] .* forceScale
        pts = [Point3f(p[1], p[2], p[3]) for p in locations]
        arrows3d!(ax, pts, dirs, color=:red)
    end

    return fig
end

end