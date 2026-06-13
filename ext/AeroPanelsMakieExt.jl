module AeroPanelsMakieExt

using GLMakie
using AeroPanels
using GeometryBasics
using StaticArrays
using LinearAlgebra

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

function add_cad_triad!(fig, main_ax)
    # 1. Create a miniature LScene superimposed in the bottom-left corner
    triad_ax = LScene(fig[1, 1]; 
        width=120, height=120, 
        halign=:left, valign=:bottom, 
        tellwidth=false, tellheight=false,
        show_axis=false
    )

    # 2. Draw the X, Y, Z orientation arrows using the updated API
    origin = [Point3f(0, 0, 0)]
    
    # 3D arrows use lineradius for the stem, and tipradius/tiplength for the cone
    arrows3d!(triad_ax, origin, [Point3f(1, 0, 0)], color=:red,   shaftradius=0.03, tipradius=0.1, tiplength=0.3)
    arrows3d!(triad_ax, origin, [Point3f(0, 1, 0)], color=:green, shaftradius=0.03, tipradius=0.1, tiplength=0.3)
    arrows3d!(triad_ax, origin, [Point3f(0, 0, 1)], color=:blue,  shaftradius=0.03, tipradius=0.1, tiplength=0.3)

    # Optional: Add axis labels
    text!(triad_ax, [Point3f(1.3, 0, 0), Point3f(0, 1.3, 0), Point3f(0, 0, 1.3)], 
          text=["X", "Y", "Z"], 
          color=[:red, :green, :blue], 
          align=(:center, :center), 
          fontsize=14)

    # 3. Synchronize camera rotation but strictly lock the zoom distance
    on(main_ax.scene.camera.view) do _
        cam_main = cameracontrols(main_ax.scene)
        cam_triad = cameracontrols(triad_ax.scene)
        
        # Get the viewing direction of the main camera
        dir = cam_main.lookat[] .- cam_main.eyeposition[]
        dir_norm = normalize(dir)
        
        # Apply the same view angle to the triad, locking the viewing distance to 4.5 units
        cam_triad.lookat[] = Vec3f(0, 0, 0)
        cam_triad.eyeposition[] = Vec3f(0, 0, 0) .- dir_norm .* 4.5
        cam_triad.upvector[] = cam_main.upvector[]
        
        update_cam!(triad_ax.scene, cam_triad)
    end
end

function AeroPanels.PlotModel(model::AeroModel; 
    plotWake=false, 
    Γp=nothing, 
    Γw=nothing, 
    plotForces=false, 
    Fa=nothing,
    kwargsLScene = (;),
    kwargsWireFrame = (;))

    fig = Figure()
    ax = LScene(fig[1, 1]; kwargsLScene...)

    add_cad_triad!(fig, ax)

    if !isnothing(Γp)
        flat_mesh_p, v_colors_p = disconnected_mesh(model.mesh, Γp)
        mesh!(ax, flat_mesh_p, color=v_colors_p, colormap=:viridis)
    else
        # Extract total panel count to preallocate the color array
        num_panels = length(faces(model.mesh))
        face_colors = fill((:lightblue, 0.75), num_panels)
        
        # Override colors for control surface panels using the stored indices
        for cs in model.controlSurfaces
            for p_idx in cs.panelIndices
                face_colors[p_idx] = (:gray, 0.75)
            end
        end
        
        flat_mesh_p, v_colors_p = disconnected_mesh(model.mesh, face_colors)
        mesh!(ax, flat_mesh_p, color=v_colors_p)
    end

    # Enforce wireframe globally to outline the transparent panels
    wireframe!(ax, model.mesh, color=:black, linewidth=1.0; kwargsWireFrame...)

    if plotWake
        if !isnothing(Γw)
            flat_mesh_w, v_colors_w = disconnected_mesh(model.wakeMesh, Γw)
            mesh!(ax, flat_mesh_w, color=v_colors_w, colormap=:viridis)
            wireframe!(ax, model.wakeMesh, color=(:black, 0.5), linewidth=0.5)
        else
            wireframe!(ax, model.wakeMesh, color=:gray, linewidth=0.5)
        end
    end

    if plotForces && !isnothing(Fa)
        c = model.modelProperties.c
        forceScale = c / maximum(norm.(Fa)) / 2
        locations = AeroPanels.AerodynamicLoadLocation(model)
        forces = Fa
        dirs = [Point3f(f[1], f[2], f[3]) for f in forces] .* forceScale
        pts = [Point3f(p[1], p[2], p[3]) for p in locations]
        arrows3d!(ax, pts, dirs, color=:red)
    end

    return fig
end

end