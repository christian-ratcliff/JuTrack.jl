using YAML
using LinearAlgebra
using Statistics
using Random
using Distributions

function load_multicharge_lattice(yaml_file::String; charges::Vector{Int}=[49, 50, 51])
    yaml_data = YAML.load_file(yaml_file)
    lattice_data = yaml_data["lattice"]
    settings = lattice_data["settings"]
    elements_data = lattice_data["elements"]
    beams_data = lattice_data["beams"]
    lines_data = lattice_data["lines"]
    main_line = lattice_data["main_line"]
    
    # Extract brho values
    brho_values = Dict{Int, Float64}()
    for (k, v) in settings["brho_values"]
        brho_values[parse(Int, string(k))] = Float64(v)
    end
    reference_charge = settings["reference_charge"]
    
    # Create ONE reference lattice (no charge-dependent scaling for SBEND)
    element_map = Dict{String, Dict{String, Any}}()
    for el in elements_data
        element_map[el["name"]] = el
    end
    
    jutrack_elements = Dict{String, Any}()
    for el in elements_data
        name = el["name"]
        type = el["type"]
        params = el["parameters"]
        
        # Create reference lattice elements (using reference charge for scaling)
        jutrack_element = create_jutrack_element(type, name, params, reference_charge, 
                                               brho_values[reference_charge], 
                                               brho_values[reference_charge], reference_charge)
        jutrack_elements[name] = jutrack_element
    end
    
    # Resolve lattice
    lattice = []
    if main_line !== nothing
        resolve_line!(lattice, main_line, lines_data, jutrack_elements)
    end
    
    # Create beams for all requested charges
    beams = Dict{Int, Beam}()
    
    # Get transformation parameters
    rest_energy = settings["IonEs"]
    kinetic_energy = settings["IonEk"]
    mass_number = settings["mass_number"]
    rf_freq = settings["rf_freq"]
    
    gamma = 1.0 + kinetic_energy / rest_energy
    beta = sqrt(1.0 - 1.0 / gamma^2)
    rf_k = 2 * pi * rf_freq / 299792458.0
    scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(rest_energy+kinetic_energy)]
    
    for charge in charges
        # Find beam data for this charge
        beam_data = nothing
        for bd in beams_data
            if bd["charge"] == charge
                beam_data = bd
                break
            end
        end
        
        if beam_data === nothing
            @warn "No beam data found for charge $charge, skipping"
            continue
        end
        
        # Create beam
        matrix_flat = beam_data["matrix"]
        matrix = reshape(matrix_flat, 6, 6)
        centroid = beam_data["centroid"]
        
        matrix_transformed = diagm(scaling) * ((matrix + matrix') ./ 2) * diagm(scaling)'
        beam = create_beam(matrix_transformed, centroid, charge, mass_number, rest_energy, kinetic_energy, brho_values, reference_charge)
        
        beams[charge] = beam
        # println("Created beam for charge $charge")
    end
    
    # Create MultiChargeBeam
    multibeam = MultiChargeBeam(beams, brho_values, reference_charge)
    
    return lattice, multibeam
end

"""
    load_lattice(yaml_file::String, charge::Int=50)

Load a JuTrack lattice from a YAML file and return the lattice elements and beam.
If `charge` is provided, it will use the corresponding beam and k1 values
for elements that depend on the charge state.
"""
function load_lattice(yaml_file::String, charge::Int=50)
    # Read the YAML file
    yaml_data = YAML.load_file(yaml_file)
    
    # Extract lattice data
    lattice_data = yaml_data["lattice"]
    settings = lattice_data["settings"]
    elements_data = lattice_data["elements"]
    beams_data = lattice_data["beams"]
    lines_data = lattice_data["lines"]
    main_line = lattice_data["main_line"]
    
    # Get properties
    mass_number = settings["mass_number"]
    rest_energy = settings["IonEs"]
    kinetic_energy = settings["IonEk"]
    rf_freq = settings["rf_freq"]
    brho_values = settings["brho_values"]
    reference_brho = settings["reference_brho"]
    reference_charge = settings["reference_charge"]

    
    
    # Check if the charge is in the available brho values
    if !haskey(brho_values, charge)
        error("Charge state $charge not found in brho_values. Available charge states: $(keys(brho_values))")
    end
    
    # Access brho directly with the charge key
    brho = brho_values[charge]
    
    # Create a mapping of element names to their definitions
    element_map = Dict{String, Dict{String, Any}}()
    for el in elements_data
        element_map[el["name"]] = el
    end
    
    # Create JuTrack elements
    jutrack_elements = Dict{String, Any}()
    
    for el in elements_data
        name = el["name"]
        type = el["type"]
        params = el["parameters"]
        
        # Create element based on type
        jutrack_element = create_jutrack_element(type, name, params, charge, brho, reference_brho, reference_charge)
        jutrack_elements[name] = jutrack_element
    end
    
    # Resolve the main line
    lattice = []
    if main_line !== nothing
        resolve_line!(lattice, main_line, lines_data, jutrack_elements)
    end
    
    # Calculate RF and beam parameters
    gamma = 1.0 + kinetic_energy / rest_energy  # Lorentz factor
    beta = sqrt(1.0 - 1.0 / (gamma * gamma))  # Velocity factor
    rf_k = 2 * pi * rf_freq / 299792458.0  # [1/m]
    
    # Define scaling factors for coordinate transformation
    scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(rest_energy+kinetic_energy)]
    
    # Find the beam for the specified charge
    beam = nothing
    for beam_data in beams_data
        if beam_data["charge"] == charge
            # Get beam data
            matrix_flat = beam_data["matrix"]
            matrix = reshape(matrix_flat, 6, 6)
            centroid = beam_data["centroid"]
            
            # Transform FLAME matrix to JuTrack matrix
            matrix_transformed = diagm(scaling) * ((matrix + matrix') ./ 2) * diagm(scaling)'
            
            # Create a beam object exactly as in the original code
            beam = create_beam(matrix_transformed, centroid, charge, mass_number, rest_energy, kinetic_energy, brho_values, reference_charge)
            break
        end
    end
    
    if beam === nothing && !isempty(beams_data)
        # If exact charge not found, use the first beam as default
        beam_data = beams_data[1]
        println("Warning: Beam with charge $charge not found, using charge $(beam_data["charge"]) as default")
        
        matrix_flat = beam_data["matrix"]
        matrix = reshape(matrix_flat, 6, 6)
        centroid = beam_data["centroid"]
        
        # Transform FLAME matrix to JuTrack matrix
        matrix_transformed = diagm(scaling) * ((matrix + matrix') ./ 2) * diagm(scaling)'
        
        beam = create_beam(matrix_transformed, centroid, beam_data["charge"], mass_number, rest_energy, kinetic_energy, brho_values, reference_charge)
    end
    
    return lattice, beam
end

function resolve_line!(lattice, line_name, lines_data, elements)
    if !haskey(lines_data, line_name)
        # Check if it's a direct element reference
        if haskey(elements, line_name)
            push!(lattice, elements[line_name])
        else
            # Instead of a dictionary, create a proper MARKER object
            # println("Warning: Element '$line_name' not found - creating a MARKER instead")
            marker = MARKER(name=line_name)
            push!(lattice, marker)
        end
        return
    end
    
    line_elements = lines_data[line_name]
    
    for item in line_elements
        element_name = item["element"]
        count = item["count"]
        
        if haskey(lines_data, element_name)
            # It's a nested line, recursively resolve it
            for _ in 1:count
                resolve_line!(lattice, element_name, lines_data, elements)
            end
        elseif haskey(elements, element_name)
            # It's a direct element
            for _ in 1:count
                push!(lattice, elements[element_name])
            end
        else
            # Create a proper MARKER object
            println("Warning: Element '$element_name' not found - creating a MARKER instead")
            marker = MARKER(name=element_name)
            # Add it the specified number of times
            for _ in 1:count
                push!(lattice, marker)
            end
        end
    end
end

"""
    create_jutrack_element(type, name, params, charge, brho)

Create a proper JuTrack element object using the appropriate element class
(DRIFT, KQUAD, SBEND, etc.) with the given parameters.
"""
function create_jutrack_element(type, name, params, charge, brho, reference_brho, reference_charge)
    brho_ratio = reference_brho / brho
    charge_ratio = charge / reference_charge
    if type == "DRIFT"
        len = get(params, "len", 0.0)
        return DRIFT(name=name, len=len)
        
    elseif type == "KQUAD"
        len = get(params, "len", 0.0)
        # Determine k1 value based on charge state
        k1 = 0.0
        if haskey(params, "k1")
            k1 = params["k1"] 
        end
        return KQUAD(name=name, len=len, k1=k1)
        
    elseif type == "SBEND"
        len = get(params, "len", 0.0)
        angle = 0.0
        if haskey(params, "angle")
            angle = params["angle"]
        elseif haskey(params, "angle_deg")
            angle = params["angle_deg"] * π / 180.0 
        end

        # Edge angles
        e1 = get(params, "e1", 0.0)
        e2 = get(params, "e2", 0.0) 
        return SBEND(name=name, len=len, angle=angle, e1=e1, e2=e2)
        
    elseif type == "ORBTRIM"
        # Handle orbit trims
        tm_xkick = get(params, "tm_xkick", 0.0)
        tm_ykick = get(params, "tm_ykick", 0.0) 
        return ORBTRIM(name=name, realpara=false, theta_x=tm_xkick, theta_y=tm_ykick)
    elseif type == "KSEXT"
        len = get(params, "len", 0.0)
        k2 = get(params, "B3", 0.0)
        return KSEXT(name=name, len=len, k2=k2)
        
    elseif type == "RFCA"
        len = get(params, "len", 0.0)
        volt = get(params, "V", 0.0)
        freq = get(params, "f", 0.0)
        lag = get(params, "lag", 0.0)
        return RFCA(name=name, len=len, volt=volt, freq=freq, lag=lag)
        
    else
        # If unrecognized element, but it has a length, create a DRIFT element
        if haskey(params, "len")
            len = get(params, "len", 0.0)
            @warn("Element type '$type' not recognized. Defaulting to DRIFT with length $len.")
            return DRIFT(name=name, len=len)
        end
        # Default to MARKER for any unrecognized element type
        return MARKER(name=name)
    end
end

"""
    create_beam(matrix, centroid, charge, mass_number, rest_energy, kinetic_energy)

Create a beam object using the JuTrack.Beam type, following the exact method from the original code.
"""
function create_beam(matrix, centroid, charge, mass_number, rest_energy, kinetic_energy, brho_values, reference_charge)
    # Calculate total energy
    total_kinetic_energy = kinetic_energy * mass_number
    total_rest_energy = rest_energy * mass_number

    # Calculate momentum deviation for this charge state
    reference_brho = brho_values[reference_charge]
    beam_brho = brho_values[Int(charge)]
    
    # Momentum deviation: δp/p = (Bρ - Bρ_ref) / Bρ_ref
    momentum_deviation = (beam_brho - reference_brho) / reference_brho
    
    Random.seed!(42 + Int(charge))  # Ensure reproducibility based on charge
    
    # Generate random particles
    nparticles = 10000
    dis = Matrix{Float64}(undef, nparticles, 6)
    for d in 1:6
        dis[:, d] .= randn(nparticles)
    end
    
    # Calculate the initial 2nd moment matrix
    moment2nd = zeros(6, 6)
    for d1 in 1:6
        for d2 in 1:6
            moment2nd[d1, d2] = mean(dis[:, d1] .* dis[:, d2])
        end
    end
    
    # Get eigendecomposition of initial matrix
    lam, u = eigen(moment2nd)
    transformation = u * diagm(1.0 ./ sqrt.(lam)) * u'
    dis = dis * transformation'
    
    # Get eigendecomposition of target matrix
    lam, u = eigen(matrix)
    
    # Create transformation matrix from eigendecomposition
    transformation = u * Diagonal(sqrt.(lam)) * u'
    
    # Apply transformation to get desired covariance
    dis = dis * transformation'

    
    
    # Calculate the final 2nd moment matrix to verify transformation
    # moment2nd = zeros(6, 6)
    for d1 in 1:6
        for d2 in 1:6
            moment2nd[d1, d2] = mean(dis[:, d1] .* dis[:, d2])
        end
    end
    
    # Calculate and display the ratio to verify accuracy
    # ratio = ifelse.((moment2nd .<= 1e-20) .& (matrix .== 0.0), 1.0, moment2nd ./ matrix)
    # println("Beam with charge $charge second moment ratio:")
    # display(ratio)
    
    # Add centroid
    dis .+= reshape(centroid, 1, 6)

    # Add energy offset for multicharge beam
    dis[:, 6] .+= momentum_deviation
    
    # Create beam with the proper JuTrack.Beam constructor
    beam = Beam(
        r=dis,  
        np=nparticles,
        energy=total_kinetic_energy,
        charge=Float64(charge),
        mass=total_rest_energy
    )
    
    # Calculate beam properties
    get_centroid!(beam)
    get_emittance!(beam)
    
    return beam
end

# If run directly, print usage
if abspath(PROGRAM_FILE) == @__FILE__
    println("Usage: include this file and call load_lattice(filename, charge)")
    println("Example: lattice, beam = load_lattice(\"jutrack_lattice.yaml\", 50)")
end
