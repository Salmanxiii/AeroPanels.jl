module FMI

using CodeTracking
using XML
using Pkg
using JuliaC
using LinearAlgebra
using Serialization
using AeroPanels
import AeroPanels: BuildFMU

include("XML.jl")
using .FMI_XML

# 1. Prebuilt Model Instance Mode
function BuildFMU(model::AbstractFMUModel, output_dir::String; 
                  fmu_name="AeroPanelsFMU", 
                  version="1.0.0",
                  extraLibraries=String[],
                  clean=true)
    println("--- AeroPanels: Starting FMI 2.0 Build (Prebuilt Model Mode) for $fmu_name ---")

    # 1. Dynamic Package Generation
    staging_dir = abspath(joinpath(output_dir, "fmu_build_staging_$(fmu_name)"))
    if isdir(staging_dir); rm(staging_dir, recursive=true, force=true); end
    mkpath(staging_dir)

    cd(staging_dir) do
        Pkg.generate(fmu_name)
    end
    
    package_dir = joinpath(staging_dir, fmu_name)
    resources_dir = joinpath(package_dir, "resources")
    mkpath(resources_dir)

    # Serialize prebuilt model instance directly to disk
    model_jls_path = joinpath(resources_dir, "model.jls")
    serialize(model_jls_path, model)
    println("Prebuilt model serialized to $(basename(model_jls_path))")

    aeropanels_path = abspath(joinpath(@__DIR__, "..", ".."))
    package_dir_esc = replace(package_dir, "\\" => "\\\\")
    
    # Setup dependencies
    defaultLibraries = ["LinearAlgebra", "StaticArrays", "SparseArrays", "GeometryBasics", "OrdinaryDiffEqTsit5", "ForwardDiff", "Serialization"]
    all_libs = unique([defaultLibraries..., extraLibraries...])
    
    setup_script = joinpath(staging_dir, "setup_project.jl")
    write(setup_script, """
    using Pkg
    using TOML
    Pkg.activate("$package_dir_esc")
    
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

    # 2. Source Generation from Template
    assets_dir = abspath(joinpath(aeropanels_path, "assets", "FMI"))
    template_path = joinpath(assets_dir, "julia", "template.jl")
    template_content = read(template_path, String)

    code_str = "# Prebuilt Model Serialization Mode"
    layout_str = """
    function create_model(resource_dir="")
        if !isempty(resource_dir)
            res_file = joinpath(resource_dir, "model.jls")
            if isfile(res_file)
                return Serialization.deserialize(res_file)
            end
        end
        candidates = [
            normpath(joinpath(@__DIR__, "..", "..", "resources", "model.jls")),
            normpath(joinpath(@__DIR__, "..", "resources", "model.jls")),
            normpath(joinpath(@__DIR__, "resources", "model.jls")),
            normpath(joinpath(pwd(), "resources", "model.jls"))
        ]
        for path in candidates
            if isfile(path)
                return Serialization.deserialize(path)
            end
        end
        error("Could not find prebuilt model serialization file model.jls in resource_dir '\$resource_dir'")
    end
    """

    content = replace(template_content, "{{FMU_NAME}}" => fmu_name)
    content = replace(content, "{{USER_CODE}}" => code_str)
    content = replace(content, "{{BUILDER_NAME}}" => "create_model")
    content = replace(content, "{{FMU_LAYOUT_DEF}}" => layout_str)
    
    write(joinpath(package_dir, "src", "$fmu_name.jl"), content)
    
    julia_assets_dir = joinpath(assets_dir, "julia")
    for f in readdir(julia_assets_dir)
        if endswith(f, ".jl") && f != "template.jl"
            cp(joinpath(julia_assets_dir, f), joinpath(package_dir, "src", f))
        end
    end

    # 3. Compilation using JuliaC Library API
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

    # 4. C-Shim Compilation
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

    # 5. Packaging
    println("Generating modelDescription.xml...")
    xml_doc = generate_xml_node(model; modelName=fmu_name, version=version)
    XML.write(joinpath(staging_dir, "modelDescription.xml"), xml_doc)

    fmu_tmp = joinpath(staging_dir, "fmu_tmp")
    bin_dest = joinpath(fmu_tmp, "binaries", "win64")
    fmu_resources_dest = joinpath(fmu_tmp, "resources")
    mkpath(bin_dest)
    mkpath(fmu_resources_dest)

    cp(joinpath(staging_dir, "modelDescription.xml"), joinpath(fmu_tmp, "modelDescription.xml"))
    cp(model_jls_path, joinpath(fmu_resources_dest, "model.jls"))

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

# 2. Parametric Factory Function Mode
function BuildFMU(builder_func::Function, output_dir::String; 
                  fmu_name="AeroPanelsFMU", 
                  version="1.0.0",
                  extraLibraries=String[],
                  clean=true)
    println("--- AeroPanels: Starting FMI 2.0 Build (Parametric Factory Mode) for $fmu_name ---")

    # Evaluate model once to check compatibility and generate XML metadata
    model = builder_func()
    if !(model isa AbstractFMUModel); error("Builder must return an AbstractFMUModel"); end

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

    staging_dir = abspath(joinpath(output_dir, "fmu_build_staging_$(fmu_name)"))
    if isdir(staging_dir); rm(staging_dir, recursive=true, force=true); end
    mkpath(staging_dir)

    cd(staging_dir) do
        Pkg.generate(fmu_name)
    end
    
    package_dir = joinpath(staging_dir, fmu_name)
    aeropanels_path = abspath(joinpath(@__DIR__, "..", ".."))
    package_dir_esc = replace(package_dir, "\\" => "\\\\")
    
    defaultLibraries = ["LinearAlgebra", "StaticArrays", "SparseArrays", "GeometryBasics", "OrdinaryDiffEqTsit5", "ForwardDiff", "Serialization"]
    all_libs = unique([defaultLibraries..., extraLibraries...])
    
    setup_script = joinpath(staging_dir, "setup_project.jl")
    write(setup_script, """
    using Pkg
    using TOML
    Pkg.activate("$package_dir_esc")
    
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

    assets_dir = abspath(joinpath(aeropanels_path, "assets", "FMI"))
    template_path = joinpath(assets_dir, "julia", "template.jl")
    template_content = read(template_path, String)

    layout_str = """
    create_model(resource_dir="") = $func_name()
    """

    content = replace(template_content, "{{FMU_NAME}}" => fmu_name)
    content = replace(content, "{{USER_CODE}}" => code_str)
    content = replace(content, "{{BUILDER_NAME}}" => func_name)
    content = replace(content, "{{FMU_LAYOUT_DEF}}" => layout_str)
    
    write(joinpath(package_dir, "src", "$fmu_name.jl"), content)
    
    julia_assets_dir = joinpath(assets_dir, "julia")
    for f in readdir(julia_assets_dir)
        if endswith(f, ".jl") && f != "template.jl"
            cp(joinpath(julia_assets_dir, f), joinpath(package_dir, "src", f))
        end
    end

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

    println("Generating modelDescription.xml...")
    xml_doc = generate_xml_node(model; modelName=fmu_name, version=version)
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
