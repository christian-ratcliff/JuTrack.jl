#use this to generate the lattice file
# julia FLAME_to_JuTrack.jl BDS_124Xe_3cs_short_v1.lat jutrack_lattice.jl

using Dates
using LinearAlgebra
using Statistics
const speed_of_light = 299792458.0  # m/s

"""
    parse_flame_lattice(filename::String, output_filename::String)

Parse a FLAME lattice file and convert it to a JuTrack lattice file.
Handles coordinate system transformation from FLAME [mm, rad, mm, rad, rad, MeV/u]
to JuTrack [m, px, m, py, m, dp].
"""
function parse_flame_lattice(filename::String, output_filename::String)
    # Read the entire file content
    file_content = read(filename, String)
    
    # Initialize collections for elements and global parameters
    element_defs = Dict{String, Dict{String, Any}}()
    global_params = Dict{String, Any}()
    
    # Extract mass number from IonChargeStates
    mass_number = 0.0  # Will be set based on file content
    charge_states_pattern = r"IonChargeStates\s*=\s*\[(.*?)\];"
    charge_match = match(charge_states_pattern, file_content)
    
    if charge_match !== nothing
        values = split(charge_match.captures[1], ",")
        for val_str in values
            val_str = strip(val_str)
            
            # Handle division expressions (like 50./124.)
            if occursin("/", val_str)
                parts = split(val_str, "/")
                if length(parts) == 2
                    try
                        # Extract the denominator which is the mass number
                        denom_str = replace(strip(parts[2]), r"\.$" => "")
                        mass_number = parse(Float64, denom_str)
                        println("Extracted mass number: $mass_number")
                        break  # We only need to extract it once
                    catch e
                        @warn "Could not extract mass number from: $val_str"
                    end
                end
            end
        end
    end
    
    if mass_number == 0.0
        # If we couldn't extract from IonChargeStates, look for it in the file
        mass_pattern = r"mass_number\s*=\s*([0-9\.]+);"
        mass_match = match(mass_pattern, file_content)
        if mass_match !== nothing
            mass_number = parse(Float64, mass_match.captures[1])
            println("Extracted mass number from explicit declaration: $mass_number")
        else
            # Set a reasonable default if we can't find it
            mass_number = 124.0  # Based on the IonChargeStates from the header
            println("Using default mass number: $mass_number")
        end
    end
    
    # Extract energy values
    rest_energy_pattern = r"IonEs\s*=\s*([0-9\.e\+\-]+);"
    rest_energy_match = match(rest_energy_pattern, file_content)
    rest_energy_per_nucleon = 931494320.0  # Default value in eV/u
    
    if rest_energy_match !== nothing
        try
            rest_energy_per_nucleon = parse(Float64, rest_energy_match.captures[1])
            println("Extracted rest energy: $rest_energy_per_nucleon eV/u")
        catch
            println("Warning: Could not parse rest energy from file, using default")
        end
    end
    
    kinetic_energy_pattern = r"IonEk\s*=\s*([0-9\.e\+\-]+);"
    kinetic_energy_match = match(kinetic_energy_pattern, file_content)
    kinetic_energy_per_nucleon = 0.0  # Will be set based on file content
    
    if kinetic_energy_match !== nothing
        try
            kinetic_energy_per_nucleon = parse(Float64, kinetic_energy_match.captures[1])
            println("Extracted kinetic energy: $kinetic_energy_per_nucleon eV/u")
        catch
            println("Warning: Could not parse kinetic energy from file")
        end
    end
    
    # Calculate relativistic parameters
    gamma = 1.0 + kinetic_energy_per_nucleon / rest_energy_per_nucleon
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beta_gamma = beta * gamma
    println("Relativistic parameters: beta=$beta, gamma=$gamma, beta*gamma=$beta_gamma")
    
    # Extract RF frequency for phase-to-length conversion
    rf_freq = 80.5e6  # Default RF frequency in Hz
    freq_pattern = r"f\s*=\s*([0-9\.e\+\-]+);"
    for freq_match in eachmatch(freq_pattern, file_content)
        try
            rf_freq = parse(Float64, freq_match.captures[1])
            println("Found RF frequency: $rf_freq Hz")
            break  # Take the first match
        catch
            continue
        end
    end
    rf_wavelength = speed_of_light / rf_freq
    println("RF wavelength: $rf_wavelength m")
    
    # Create the coordinate transformation matrix
    # FLAME: [x(mm), x'(rad), y(mm), y'(rad), φ(rad), dE_k(MeV/u)]
    # JuTrack: [x(m), px, y(m), py, z(m), dp/p]
    T = zeros(Float64, 6, 6)
    scaling = zeros(Float64, 6)
    scaling = [1e-3, 1.0, 1e-3, 1.0, rf_wavelength * beta / (2.0 * π),  1.0 / (beta^2 * kinetic_energy_per_nucleon * 1.0e-6) ] 
    # scaling .= 1.0
    T = diagm(scaling)
    
    # T = scaling * scaling'

    
    println("Coordinate transformation matrix T:")
    display(T)
    println()
    
    # Extract S matrices for beam envelopes
    beam_envelopes = Dict{String, Matrix{Float64}}()
    matrix_names = ["S0", "S1", "S2"]
    
    for name in matrix_names
        println("Looking for matrix $name")
        pattern = Regex("$name\\s*=\\s*\\[(.*?)\\];", "s")
        m = match(pattern, file_content)
        
        if m !== nothing
            println("Found matrix $name")
            matrix_data = m.captures[1]
            
            # Create a 6x6 matrix
            matrix = zeros(Float64, 6, 6)
            
            # Split into rows
            rows = split(matrix_data, "\n")
            rows = [strip(row) for row in rows if !isempty(strip(row))]
            
            for i in 1:min(6, length(rows))
                # Split row into values
                values = split(rows[i], ",")
                values = [strip(v) for v in values if !isempty(strip(v))]
                
                for j in 1:min(6, length(values))
                    try
                        # Handle scientific notation and other numeric formats
                        val_str = values[j]
                        val = nothing
                        try
                            val = parse(Float64, val_str)
                        catch
                            try
                                val = eval(Meta.parse(val_str))
                            catch
                                println("Warning: Could not parse value '$val_str' at position [$i,$j]")
                                val = 0.0
                            end
                        end
                        matrix[i, j] = val
                    catch e
                        println("Error parsing value at position [$i,$j]: $e")
                        matrix[i, j] = 0.0
                    end
                end
            end
            # display(cor(matrix))
            
            # Apply coordinate transformation: Σ' = T·Σ·Tᵀ
            jutrack_matrix = T * matrix * T'
            # Store the transformed matrix
            beam_envelopes[name] = jutrack_matrix
            
            println("Transformed matrix $name to JuTrack coordinates")
        else
            println("Warning: Could not find matrix $name")
        end
    end
    
    # Extract and transform BaryCenter vectors
    barycenters = Dict{String, Vector{Float64}}()
    
    for i in 0:2
        name = "BaryCenter$i"
        pattern = Regex("$name\\s*=\\s*\\[(.*?)\\];", "s")
        m = match(pattern, file_content)
        if m !== nothing
            # Extract the values as a string
            values_str = m.captures[1]

            # Split the values string by commas and strip spaces
            values_str_split = split(values_str, ",")
            
            # Evaluate each value individually and convert to a number
            evaluated_values = Float64[]
            for value in values_str_split
                value = strip(value)  # Remove any extra spaces
                try
                    # Evaluate the value if it's a mathematical expression
                    eval_value = eval(Meta.parse(value))
                    push!(evaluated_values, eval_value)
                catch e
                    println("Error evaluating value: $value")
                    push!(evaluated_values, 0.0)  # Use default if parsing fails
                end
            end
            
            # Remove the 7th element if it exists (the weight/fraction)
            if length(evaluated_values) >= 7
                deleteat!(evaluated_values, 7)
            end
            
            # Pad to length 6 if shorter
            while length(evaluated_values) < 6
                push!(evaluated_values, 0.0)
            end
            
            # Convert FLAME centroid to JuTrack centroid using the transformation matrix
            flame_centroid = evaluated_values[1:6]
            jutrack_centroid = scaling .* flame_centroid
            
            # Store the transformed centroid
            barycenters[name] = jutrack_centroid
            
            println("Transformed $name to JuTrack coordinates: $jutrack_centroid")
        else
            println("Warning: Could not find $name")
        end
    end
    
    # Extract charge states and NCharge values
    charge_states = Float64[]
    if charge_match !== nothing
        println("Found IonChargeStates: $(charge_match.captures[1])")
        values = split(charge_match.captures[1], ",")
        for val_str in values
            val_str = strip(val_str)
            
            # Handle division expressions (like 50./124.)
            if occursin("/", val_str)
                parts = split(val_str, "/")
                if length(parts) == 2
                    try
                        # Parse numerator as a float
                        numer = parse(Float64, replace(strip(parts[1]), r"\.$" => ""))
                        charge_states = push!(charge_states, numer)
                    catch e
                        @warn "Could not parse charge state fraction: $val_str"
                    end
                end
            else
                # Try direct parsing
                try
                    push!(charge_states, parse(Float64, val_str))
                catch e
                    @warn "Could not parse charge state: $val_str"
                end
            end
        end
        println("Parsed charge states: $charge_states")
    else
        println("Warning: Could not find IonChargeStates")
    end
    
    # Extract NCharge values
    ncharge = Float64[]
    ncharge_pattern = r"NCharge\s*=\s*\[(.*?)\];"
    ncharge_match = match(ncharge_pattern, file_content)
    if ncharge_match !== nothing
        println("Found NCharge: $(ncharge_match.captures[1])")
        values = split(ncharge_match.captures[1], ",")
        for val_str in values
            val_str = strip(val_str)
            try
                push!(ncharge, parse(Float64, val_str))
            catch e
                @warn "Could not parse NCharge value: $val_str"
            end
        end
        println("Parsed NCharge values: $ncharge")
    else
        println("Warning: Could not find NCharge")
    end
    
    # Extract element definitions and maintain order
    element_order = String[]
    element_pattern = r"([^;:]+(?::[^;:]+)?)\s*:\s*([^,;]+)(?:,\s*(.*?))?;"
    for m in eachmatch(element_pattern, file_content)
        element_name = String(strip(m[1]))
        element_type = String(strip(m[2]))
        
        # Skip LINE definitions, comments, source, or explanatory text
        if occursin("LINE", element_type) || 
           element_type == "source" || 
           startswith(element_name, "#") || 
           occursin("Note", element_name) ||
           occursin("Units are", element_name) ||
           occursin("is:", element_name) ||
           occursin("are:", element_name) ||
           length(element_name) > 50  # Likely explanatory text, not an element name
            continue
        end
        
        # Add to the ordered list if not already added
        if !(element_name in element_order)
            push!(element_order, element_name)
        end
        
        # Parse parameters
        element_params = Dict{String, Any}()
        element_params["type"] = element_type
        element_params["name"] = element_name
        
        if m[3] !== nothing
            params_str = m[3]
            for param in split(params_str, ",")
                param = strip(param)
                if occursin("=", param)
                    param_parts = split(param, "=", limit=2)
                    param_name = String(strip(param_parts[1]))
                    param_value_str = String(strip(param_parts[2]))
                    
                    if param_name == "bg"  # Skip beam gamma parameter
                        continue
                    end

                    # Try to convert to appropriate type
                    param_value = try
                        parse(Float64, param_value_str)
                    catch
                        # Remove quotes and try numeric evaluation
                        clean_value = replace(param_value_str, "\"" => "")
                        try
                            if occursin(r"[0-9+\-*/\^\.e]", clean_value) && !occursin(r"[a-zA-Z]", clean_value)
                                eval(Meta.parse(clean_value))
                            else
                                clean_value
                            end
                        catch
                            clean_value
                        end
                    end
                    
                    element_params[param_name] = param_value
                end
            end
        end
        
        element_defs[element_name] = element_params
    end
    
    # Extract energy and mass values with proper unit conversion
    total_kinetic_energy = kinetic_energy_per_nucleon * mass_number  # Total kinetic energy in eV
    total_rest_energy = rest_energy_per_nucleon * mass_number  # Total rest energy in eV
    total_energy = total_kinetic_energy + total_rest_energy  # Total energy in eV
    
    # Now generate the JuTrack lattice file
    open(output_filename, "w") do file
        # Write header
        write(file, """
        # JuTrack lattice file converted from FLAME format
        # Original file: $filename
        # Generated on: $(Dates.now())
        # Mass number: $mass_number
        # Coordinate system conversion: FLAME [mm, rad, mm, rad, rad, MeV/u] -> JuTrack [m, px, m, py, m, dp]
        # Relativistic factors: beta=$beta, gamma=$gamma, beta*gamma=$beta_gamma
        # RF frequency: $rf_freq Hz, wavelength: $rf_wavelength m
        
        include("../../JuTrack.jl")
        using .JuTrack
        using LinearAlgebra
        using Random
        using Distributions
        using CairoMakie
        
        """
        )
        
        # Write function to create the lattice
        write(file, """
        function create_lattice()
            lattice = []
        
            # Define elements
        """)
        
        # Define all elements in order
        for element_name in element_order
            if haskey(element_defs, element_name)
                params = element_defs[element_name]
                
                # Generate a safe name for Julia
                safe_name = replace(element_name, r"[: #\.]" => "_")
                safe_name = replace(safe_name, r"[\(\)]" => "")
                
                # Map FLAME element type to JuTrack type
                jutrack_type = get_jutrack_type(params["type"])
                
                # Generate parameter list
                param_str = "name=\"$(element_name)\""
                
                for (param_name, param_value) in params
                    if param_name != "type" && param_name != "name"
                        # Map parameter names and apply unit conversions if needed
                        jutrack_param = get_jutrack_param(param_name, param_value, params["type"], mass_number)
                        if !isnothing(jutrack_param)
                            param_str *= ", $(jutrack_param)"
                        end
                    end
                end
                
                # Create element
                write(file, "    $safe_name = $jutrack_type($param_str)\n")
                
                # Add element to the lattice
                write(file, "    push!(lattice, $safe_name)\n")
            end
        end
        
        # Write matrices
        write(file, "\n    # Define beam envelope matrices\n")
        for (matrix_name, matrix) in beam_envelopes
            write(file, "    $(matrix_name)_matrix = [\n")
            for i in 1:6
                write(file, "        ")
                for j in 1:6
                    write(file, "$(matrix[i, j]) ")
                end
                write(file, ";\n")
            end
            write(file, "    ]\n\n")
        end
        
        # Write beam creation code
        write(file, "    # Create beam objects for each charge state\n")
        
        for i in 1:min(3, length(charge_states))
            # Get barycenter data
            bary_key = "BaryCenter$(i-1)"
            centroid = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Default centroid
            
            if haskey(barycenters, bary_key)
                centroid = barycenters[bary_key]
            end
            
            # Get charge value and weight
            charge_value = charge_states[i]
            weight_value = (i <= length(ncharge)) ? ncharge[i] : 1.0
            
            # Create beam with transformed coordinates
            write(file, """
            beam$i = create_beam_from_envelope_matrix(
                S$(i-1)_matrix, 
                centroid = [$(centroid[1]), $(centroid[2]), $(centroid[3]), $(centroid[4]), $(centroid[5]), $(centroid[6])], 
                nparticles = 10000, 
                energy = $total_kinetic_energy, 
                charge = $charge_value, 
                mass = $total_rest_energy)
            
            """)
        end
        
        write(file, "    return lattice, beam1, beam2, beam3, S0_matrix, S1_matrix, S2_matrix\n")
        write(file, "end\n\n")
        
        # Add helper functions
        write(file, """
        # Helper function to create a beam from an envelope matrix
        function create_beam_from_envelope_matrix(envelope_matrix::Matrix{Float64}; 
                                        centroid::Vector{Float64} = zeros(6),
                                        nparticles::Int = 10000, 
                                        energy::Float64 = 1e9, 
                                        charge::Float64 = -1.0, 
                                        mass::Float64 = 931494320.0, 
                                        seed::Int = 42)
            # Validate inputs
            if size(envelope_matrix) != (6, 6)
                error("Envelope matrix must be 6x6")
            end
            
            if length(centroid) != 6
                error("Centroid must be a 6-element vector")
            end
            
            Random.seed!(seed)
            
            # Ensure the envelope matrix is symmetric
            envelope_matrix = (envelope_matrix + envelope_matrix')/2

            # Generate random samples with the correct covariance
            # Method 1: Use Cholesky decomposition if matrix is positive definite
            L = try
                cholesky(envelope_matrix).L
            catch
                # Method 2: Use eigendecomposition as fallback
                eigen_decomp = eigen(envelope_matrix)
                eigen_vals = eigen_decomp.values
                eigen_vecs = eigen_decomp.vectors
                # Ensure all eigenvalues are positive (for numerical stability)
                eigen_vals = max.(eigen_vals, 1e-20)
                
                # Compute L = V * sqrt(D)
                eigen_vecs * Diagonal(sqrt.(eigen_vals))
            end
            
            # Generate standard normal random variables
            Z = randn(nparticles, 6)
            
            particles = Z * transpose(L)
            particles .+= reshape(centroid, 1, 6)
            
            # Create beam with the generated particles
            beam = Beam(r=particles, np=nparticles, energy=energy, charge=charge, mass=mass)
            get_centroid!(beam)
            get_emittance!(beam)

            
            return beam
        end
        
        
        # Run the simulation
        function run_simulation()
            lattice, beam1, beam2, beam3, S0_matrix, S1_matrix, S2_matrix = create_lattice()
            
            linepass!(lattice, beam1)
            get_emittance!(beam1)
            get_centroid!(beam1)

            linepass!(lattice, beam2)
            get_emittance!(beam2)
            get_centroid!(beam2)

            linepass!(lattice, beam3)
            get_emittance!(beam3)
            get_centroid!(beam3)

            println("Final beam1 parameters:")
            println("Centroid: ", beam1.emittance)
            println("Emittance: ", beam1.centroid)

            println("Final beam2 parameters:")
            println("Centroid: ", beam2.emittance)
            println("Emittance: ", beam2.centroid)

            println("Final beam3 parameters:")
            println("Centroid: ", beam3.emittance)
            println("Emittance: ", beam3.centroid)
            
            return lattice, beam1, beam2, beam3, S0_matrix, S1_matrix, S2_matrix
        end
        
        function visualize_beam_properties(beam::Beam)
            # Create a multi-panel plot to visualize beam properties
            
            # Phase space plots (x-px, y-py, z-dp)
            p1 = scatter(beam.r[:,1], beam.r[:,2], 
                        markersize=1, markerstrokewidth=0, alpha=0.5,
                        title="x-px phase space", xlabel="x [m]", ylabel="px")
                        
            p2 = scatter(beam.r[:,3], beam.r[:,4], 
                        markersize=1, markerstrokewidth=0, alpha=0.5,
                        title="y-py phase space", xlabel="y [m]", ylabel="py")
                        
            p3 = scatter(beam.r[:,5], beam.r[:,6], 
                        markersize=1, markerstrokewidth=0, alpha=0.5,
                        title="z-dp phase space", xlabel="z [m]", ylabel="dp/p")
            
            # Projections (histograms)
            p4 = histogram(beam.r[:,1], bins=50, alpha=0.7, title="x distribution", xlabel="x [m]")
            p5 = histogram(beam.r[:,3], bins=50, alpha=0.7, title="y distribution", xlabel="y [m]")
            p6 = histogram(beam.r[:,5], bins=50, alpha=0.7, title="z distribution", xlabel="z [m]")
            
            # Combine plots
            plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(900,600), legend=false)
        end
        
        # Create the lattice and beams
        lattice, beam1, beam2, beam3, S0_matrix, S1_matrix, S2_matrix = create_lattice();
        
        begin
            # Print validation information
            println("Validation results:")
            println("beam1 covariance ratio:")
            ratio1 = cov(beam1.r) ./ S0_matrix
            display(ratio1)
            
            println("\\nbeam2 covariance ratio:")
            ratio2 = cov(beam2.r) ./ S1_matrix
            display(ratio2)
            
            println("\\nbeam3 covariance ratio:")
            ratio3 = cov(beam3.r) ./ S2_matrix
            display(ratio3)
        end
        lattice, beam1, beam2, beam3, S0_matrix, S1_matrix, S2_matrix = run_simulation();
        # Run visualization
        p1 = visualize_beam_properties(beam1)
        p2 = visualize_beam_properties(beam2)
        p3 = visualize_beam_properties(beam3)
        

        """)
    end
    
    println("JuTrack lattice file generated: $output_filename")
    return true
end

function get_jutrack_type(flame_type::String)
    # Make sure the input is a String, not a SubString
    flame_type = String(flame_type)
    
    type_map = Dict(
        "source" => "MARKER",
        "bpm" => "MARKER",
        "drift" => "DRIFT",
        "orbtrim" => "ORBTRIM",
        "quadrupole" => "KQUAD",
        "marker" => "MARKER",
        "sbend" => "SBEND",
        "sextupole" => "KSEXT",
        "rfcavity" => "RFCA",
        "stripper" => "STRIPPER",
        "solenoid" => "SOLENOID"
    )
    
    return get(type_map, flame_type, "MARKER")  # Default to MARKER if type not recognized
end


function get_jutrack_param(param_name::String, param_value, flame_type::String, mass_number::Float64)
    # Make sure param_name is a String, not a SubString
    param_name = String(param_name)
    flame_type = String(flame_type)
    
    # Skip these parameters
    if param_name in ["aper", "vector_variable", "matrix_variable"]
        return nothing
    end
    
    # Parameter mappings
    param_map = Dict(
        "L" => "len",
        "B2" => "k1",
        "B3" => "k2",
        "B4" => "k3",
        "phi" => "angle",
        "phi1" => "e1",
        "phi2" => "e2",
        "tm_xkick" => "tm_xkick",
        "tm_ykick" => "tm_ykick",
        "realpara" => "realpara",
        "V" => "volt",
        "f" => "freq",
        "lag" => "lag",
        "phi_s" => "phis",
        "B" => "ks",
        "IonZ" => "IonZ",
        "IonMass" => "IonMass",
        "IonProton" => "IonProton"
    )
    
    # Special case for realpara
    if param_name == "realpara" && param_value == 1.0
        return "realpara=true"
    end
    
    # Special case for energy parameters - convert from eV/u to eV
    if param_name in ["IonEs", "IonEk"] && typeof(param_value) <: Number
        return "$(get(param_map, param_name, param_name))=$(param_value * mass_number)"
    end
    
    # Return mapped parameter name and value
    mapped_name = get(param_map, param_name, param_name)
    return "$(mapped_name)=$(param_value)"
end

# Process command line arguments
if length(ARGS) < 1 || length(ARGS) > 2
    println("Usage: julia FLAME_to_JuTrack.jl <input_flame_file> [output_jutrack_file]")
    println("       If output file is not specified, it will default to 'jutrack_lattice.jl'")
    exit(1)
end

input_file = ARGS[1]
output_file = length(ARGS) == 2 ? ARGS[2] : "jutrack_lattice.jl"

# Check if input file exists
if !isfile(input_file)
    println("Error: Input file '$input_file' does not exist.")
    exit(1)
end

# Run the conversion
println("Converting FLAME lattice file '$input_file' to JuTrack format...")

try
    success = parse_flame_lattice(input_file, output_file)
    if success
        println("Conversion successful! JuTrack lattice file saved to '$output_file'")
        exit(0)
    else
        println("Error: Conversion failed.")
        exit(1)
    end
catch e
    println("Error during conversion: $e")
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    exit(1)
end