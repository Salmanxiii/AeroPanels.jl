module FMI

using CodeTracking
using XML
using Pkg
using JuliaC
using LinearAlgebra
using AeroPanels
import AeroPanels: BuildFMU

include("XML.jl")
using .FMI_XML

function BuildFMU(builder_func::Function, output_dir::String; 
                  fmu_name="AeroPanelsFMU", 
                  version="1.0.0",
                  extraLibraries=String[],
                  clean=true)
    println("--- AeroPanels: Starting FMI 2.0.5 Build for $fmu_name ---")

    # 1. Validation and Source Extraction
    m = methods(builder_func)
    if length(m) == 0; error("Function has no methods"); end
    
    code = definition(m[1])
    code_str = ""
    if code === nothing
        println("Warning: CodeTracking failed, trying manual fallback...")
        meth = m[1]
        file, line = string(meth.file), meth.line
        if isfile(file)
            lines = readlines(file)
            func_lines = String[]
            for i in line:length(lines)
                push!(func_lines, lines[i])
                if i > line && (startswith(strip(lines[i]), "end")); break; end
            end
            code_str = join(func_lines, "\n")
        else
            error("Source file $file not found.")
        end
    else
        code_str = string(code)
    end
    func_name = string(nameof(builder_func))

    # Get model once to check type and dimensions
    model = builder_func()
    if !(model isa UnsteadyAeroModel2D); error("Builder must return an UnsteadyAeroModel2D"); end
    nx = NumberOfStates(model)
    nCtrl = length(model.controlSurfaces)
    nVert = model.sizes.totalVertices
    nmp = length(model.monitorPoints)

    # 2. Dynamic Package Generation
    staging_dir = abspath(joinpath(output_dir, "fmu_build_staging_$(fmu_name)"))
    if isdir(staging_dir); rm(staging_dir, recursive=true, force=true); end
    mkpath(staging_dir)

    cd(staging_dir) do
        Pkg.generate(fmu_name)
    end
    
    package_dir = joinpath(staging_dir, fmu_name)
    aeropanels_path = abspath(joinpath(@__DIR__, "..", ".."))
    package_dir_esc = replace(package_dir, "\\" => "\\\\")
    
    # Setup dependencies
    defaultLibraries = ["LinearAlgebra", "StaticArrays", "SparseArrays", "GeometryBasics", "OrdinaryDiffEqTsit5"]
    all_libs = unique([defaultLibraries..., extraLibraries...])
    
    setup_script = joinpath(staging_dir, "setup_project.jl")
    write(setup_script, """
    using Pkg
    using TOML
    Pkg.activate("$package_dir_esc")
    
    # Update Project.toml metadata manually via TOML first
    project_path = Pkg.project().path
    project_data = TOML.parsefile(project_path)
    project_data["version"] = "$version"
    project_data["authors"] = ["AeroPanels.jl"]
    open(project_path, "w") do io
        TOML.print(io, project_data)
    end
    
    Pkg.develop(path="$(replace(aeropanels_path, "\\" => "\\\\"))")
    Pkg.add([$(join(map(s -> "\"$s\"", all_libs), ", "))])
    
    Pkg.instantiate()
    """)
    run(`julia +release $setup_script`)

    # 3. Source Generation from Template
    assets_dir = abspath(joinpath(aeropanels_path, "assets", "FMI"))
    template_path = joinpath(assets_dir, "julia", "template.jl")
    template_content = read(template_path, String)
    
    # Sequential flat ValueReference offsets for FMI 2.0 scalar variables
    vr_x = 1
    vr_der_x = nx + 1
    vr_vb = 2 * nx + 1
    vr_wb = 2 * nx + 4
    vr_ab = 2 * nx + 7
    vr_dwb = 2 * nx + 10
    vr_rho = 2 * nx + 13
    vr_cs_start = 2 * nx + 14
    
    vr_rs = 2 * nx + 14 + 3 * nCtrl
    vr_vs = vr_rs + 3 * nVert
    vr_as = vr_vs + 3 * nVert
    
    vr_out_start = 2 * nx + 14 + 3 * nCtrl + (nVert > 0 ? 9 * nVert : 0)
    vr_Fb = vr_out_start
    vr_Mb = vr_out_start + 3
    vr_mp_start = vr_out_start + 6

    layout_str = """
    const FMU_LAYOUT = (
        nx = $nx,
        nCtrl = $nCtrl,
        nVert = $nVert,
        nmp = $nmp,
        
        dim_vb = 3,
        dim_wb = 3,
        dim_ab = 3,
        dim_dwb = 3,
        dim_Fb = 3,
        dim_Mb = 3,
        
        vr_x = $vr_x,
        vr_der_x = $vr_der_x,
        vr_vb = $vr_vb,
        vr_wb = $vr_wb,
        vr_ab = $vr_ab,
        vr_dwb = $vr_dwb,
        vr_rho = $vr_rho,
        
        vr_cs_start = $vr_cs_start,
        vr_rs = $vr_rs,
        vr_vs = $vr_vs,
        vr_as = $vr_as,
        
        vr_Fb = $vr_Fb,
        vr_Mb = $vr_Mb,
        vr_mp_start = $vr_mp_start,
    )
    """

    content = replace(template_content, "{{FMU_NAME}}" => fmu_name)
    content = replace(content, "{{USER_CODE}}" => code_str)
    content = replace(content, "{{BUILDER_NAME}}" => func_name)
    content = replace(content, "{{FMU_LAYOUT_DEF}}" => layout_str)
    
    # Write main module file
    write(joinpath(package_dir, "src", "$fmu_name.jl"), content)
    
    # Copy FMI Julia headers from assets/FMI/julia/ into src/
    julia_assets_dir = joinpath(assets_dir, "julia")
    for f in readdir(julia_assets_dir)
        if endswith(f, ".jl") && f != "template.jl"
            cp(joinpath(julia_assets_dir, f), joinpath(package_dir, "src", f))
        end
    end

    # 4. Compilation using JuliaC Library API
    build_dir = joinpath(staging_dir, "build")
    
    img = ImageRecipe(
        output_type = "--output-lib",
        file        = package_dir,
        trim_mode   = "no",
        add_ccallables = true,
        verbose     = true,
    )

    link = LinkRecipe(
        image_recipe = img,
        outname      = joinpath(build_dir, "bin", "fmu_engine"),
        rpath        = "@bundle",
    )

    bun = BundleRecipe(
        link_recipe = link,
        output_dir  = build_dir, 
    )

    compile_products(img)
    link_products(link)
    bundle_products(bun)

    # 5. C-Shim Compilation
    c_assets_dir = joinpath(assets_dir, "c")
    shim_source = joinpath(c_assets_dir, "fmi2_shim.c")
    shim_dest = joinpath(build_dir, "bin", fmu_name * ".dll")

    mingw_bin = joinpath(homedir(), ".julia", "artifacts", "b17bda08a19173572926f43a48aad5ef3d845e7c", "extracted_files", "mingw64", "bin")
    old_path = get(ENV, "PATH", "")
    ENV["PATH"] = mingw_bin * ";" * old_path
    try
        run(`gcc -shared -o $shim_dest $shim_source -I$c_assets_dir`)
    finally
        ENV["PATH"] = old_path
    end

    # 6. Packaging
    println("Generating modelDescription.xml...")
    nCtrl = length(model.controlSurfaces)
    csNames = [cs.name for cs in model.controlSurfaces]
    nVert = model.sizes.totalVertices
    nmp = length(model.monitorPoints)
    mpNames = [mp.name for mp in model.monitorPoints]
    xml_doc = generate_xml_node(fmu_name, nx, nCtrl, csNames, nVert, nmp, mpNames, version)
    XML.write(joinpath(staging_dir, "modelDescription.xml"), xml_doc)

    fmu_tmp = joinpath(staging_dir, "fmu_tmp")
    bin_dest = joinpath(fmu_tmp, "binaries", "win64")
    mkpath(bin_dest)
    cp(joinpath(staging_dir, "modelDescription.xml"), joinpath(fmu_tmp, "modelDescription.xml"))
    for f in readdir(joinpath(build_dir, "bin"))
        cp(joinpath(build_dir, "bin", f), joinpath(bin_dest, f))
    end

    fmu_path = abspath(joinpath(output_dir, fmu_name * ".fmu"))
    if isfile(fmu_path); rm(fmu_path, force=true); end
    zip_path = fmu_path * ".zip"
    run(`powershell -NoProfile -Command "Compress-Archive -Force -Path '$fmu_tmp/*' -DestinationPath '$zip_path'; Rename-Item -Path '$zip_path' -NewName '$(fmu_name).fmu'"`)

    if clean; rm(staging_dir, recursive=true, force=true); end
    println("--- AeroPanels: FMI Build Complete! ---")
    println("FMU created at: $fmu_path")
end

end
