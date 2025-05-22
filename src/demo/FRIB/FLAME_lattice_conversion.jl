using YAML
using Dates
using LinearAlgebra
using Statistics

const speed_of_light = 299792458.0  # m/s

function zero_non_diagonal_blocks(matrix)
    # Check that we have a square matrix
    rows, cols = size(matrix)
    if rows != cols || rows % 2 != 0
        error("Matrix must be square with even dimensions")
    end
    
    # Create a copy of the matrix
    result = copy(matrix)
    
    # Size of sub-matrices
    block_size = 2
    
    # Zero out all elements not in diagonal blocks
    for i in 1:rows
        for j in 1:cols
            # Calculate which block this element belongs to
            block_i = ceil(Int, i/block_size)
            block_j = ceil(Int, j/block_size)
            
            # If not in a diagonal block, set to zero
            if block_i != block_j
                result[i, j] = 0
            end
        end
    end
    
    return result
end

function parse_flame_lattice(filename::String, output_filename::String)
    # Read the entire file content
    file_content = read(filename, String)
    
    # Initialize collections for elements and global parameters
    element_defs = Dict{String, Dict{String, Any}}()
    global_params = Dict{String, Any}()
    
    # Extract mass number from IonChargeStates
    mass_number = 0.0
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
                        break
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
        mass_number = parse(Float64, mass_match.captures[1])
        println("Extracted mass number from explicit declaration: $mass_number")

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
    kinetic_energy_per_nucleon = 227050000.0  # Default value
    
    if kinetic_energy_match !== nothing
        try
            kinetic_energy_per_nucleon = parse(Float64, kinetic_energy_match.captures[1])
            println("Extracted kinetic energy: $kinetic_energy_per_nucleon eV/u")
        catch
            println("Warning: Could not parse kinetic energy from file, using default")
        end
    end
    
    # Calculate relativistic parameters
    gamma = 1.0 + kinetic_energy_per_nucleon / rest_energy_per_nucleon
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beta_gamma = beta * gamma
    println("Relativistic factors: beta=$beta, gamma=$gamma, beta*gamma=$beta_gamma")
    
    # Define brho values for different charge states
    brho_values = Dict(
        51 => 5.57963986,
        50 => 5.69860938,
        49 => 5.81944925
    )
    
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
    rf_k = 2 * pi * rf_freq / speed_of_light  # [1/m]
    println("RF wavelength: $rf_wavelength m")
    
    # Create the scaling array for coordinate transformation
    # FLAME: [x(mm), x'(rad), y(mm), y'(rad), φ(rad), dE_k(MeV/u)]
    # JuTrack: [x(m), px, y(m), py, z(m), dp/p]
    scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(rest_energy_per_nucleon+kinetic_energy_per_nucleon)]
    
    println("Coordinate transformation scaling:")
    display(scaling)
    println()
    
    # Extract S matrices for beam envelopes
    raw_matrices = Dict{String, Matrix{Float64}}()
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
            
            # Check eigenvalues
            lam, u = eigen(matrix)
            if any(lam .<= 0)
                println("Warning: Detected negative eigenvalues in the $name matrix, resetting matrix to zero cross terms")
                matrix = zero_non_diagonal_blocks(matrix)
            end
            
            # Store the raw matrix
            raw_matrices[name] = matrix
            println("Stored raw matrix $name")
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
            
            # Convert FLAME centroid to JuTrack centroid using the scaling
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
            
            # Handle division expressions
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
    
    # Extract element definitions
    element_pattern = r"^[ \t]*([A-Za-z0-9_:\.]+)[ \t]*:[ \t]*([A-Za-z0-9_]+)(?:,[ \t]*(.*?))?;[ \t]*$"m
    
    for m in eachmatch(element_pattern, file_content)
        element_name = String(strip(m[1]))
        element_type = String(strip(m[2]))
        
        # Skip LINE definitions, comments, source, or explanatory text
        if occursin("LINE", element_type) || 
            element_type == "source" || 
            startswith(element_name, "#") || 
            occursin("Note", element_name) ||
            occursin("#", element_name) ||
            occursin("Units are", element_name) ||
            occursin("is:", element_name) ||
            occursin("are:", element_name) ||
            length(element_name) > 50  # Likely explanatory text, not an element name
            continue
        end
        
        element_name = split(element_name, '\n')[end]  # Take the last line if multi-line
        element_name = strip(element_name)

        if isempty(element_name) || startswith(element_name, "#")
            continue
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
    
    # Extract LINE definitions to handle element repetitions
    line_defs = Dict{String, Vector{Dict{String, Any}}}()
    line_pattern = r"([A-Za-z0-9_]+):\s*LINE\s*=\s*\((.*?)\);"s
    
    for m in eachmatch(line_pattern, file_content)
        line_name = strip(m[1])
        line_content = strip(m[2])
        
        # Parse line content
        line_elements = Vector{Dict{String, Any}}()
        for element_str in split(line_content, ",")
            element_str = strip(element_str)
            if isempty(element_str)
                continue
            end
            
            # Check for repetition (element*N)
            rep_match = match(r"(.*)\*\s*(\d+)", element_str)
            if rep_match !== nothing
                element_name = strip(rep_match[1])
                count = parse(Int, rep_match[2])
            else
                element_name = element_str
                count = 1
            end
            
            push!(line_elements, Dict("element" => element_name, "count" => count))
        end
        
        line_defs[line_name] = line_elements
    end
    
    # Find the main line that is used (USE: line_name)
    use_pattern = r"USE:\s*([A-Za-z0-9_]+);"
    use_match = match(use_pattern, file_content)
    main_line = use_match !== nothing ? strip(use_match[1]) : nothing

    reference_charge = 50
    reference_brho = brho_values[reference_charge]
    
    # Create the YAML structure
    yaml_dict = Dict{String, Any}(
        "lattice" => Dict{String, Any}(
            "settings" => Dict{String, Any}(
                "mass_number" => mass_number,
                "IonEs" => rest_energy_per_nucleon,
                "IonEk" => kinetic_energy_per_nucleon,
                "rf_freq" => rf_freq,
                "brho_values" => brho_values,
                "reference_charge" => reference_charge,
                "reference_brho" => reference_brho
            ),
            "elements" => [],
            "beams" => [],
            "lines" => line_defs,
            "main_line" => main_line
        )
    )
    
    # Add all elements
    for (name, params) in element_defs
        # Map element type to JuTrack type
        jutrack_type = get_jutrack_type(params["type"])
        
        # Convert parameters
        jutrack_params = Dict{String, Any}()
        for (param_name, param_value) in params
            if param_name != "type" && param_name != "name" && param_name != "bg"
                # Handle special case for quadrupole strength
                if params["type"] == "quadrupole" && param_name == "B2"
                    # Store the raw B2 value and the k1 values for each charge state
                    jutrack_params["B2"] = param_value
                    jutrack_params["k1"] = param_value / reference_brho
                    
                # Handle angle conversion for bends
                elseif params["type"] == "sbend" && (param_name == "phi" || param_name == "angle")
                    jutrack_params["angle_deg"] = param_value  # Store original angle in degrees
                    jutrack_params["angle"] = param_value * π / 180.0  # Convert to radians
                # Handle edge angles for bends
                elseif params["type"] == "sbend" && (param_name == "phi1" || param_name == "phi2")
                    jutrack_params[param_name == "phi1" ? "e1" : "e2"] = param_value * π / 180.0
                # Standard length parameter
                elseif param_name == "L"
                    jutrack_params["len"] = param_value
                # Handle orbit corrector kicks 
                elseif params["type"] == "orbtrim" && (param_name == "tm_xkick" || param_name == "tm_ykick")
                    jutrack_params[param_name] = param_value / reference_brho # Don't scale kicks
                # All other parameters
                else
                    jutrack_params[param_name] = param_value
                end
            end
        end
        
        push!(yaml_dict["lattice"]["elements"], Dict(
            "name" => name,
            "type" => jutrack_type,
            "parameters" => jutrack_params
        ))
    end
    
    # Add beams
    for i in 1:min(3, length(charge_states))
        charge_value = charge_states[i]
        bary_key = "BaryCenter$(i-1)"
        matrix_key = "S$(i-1)"
        
        centroid = haskey(barycenters, bary_key) ? barycenters[bary_key] : [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        matrix = haskey(raw_matrices, matrix_key) ? raw_matrices[matrix_key] : zeros(6, 6)
        
        # Make matrix 2D for YAML storage
        matrix_flat = reshape(matrix, 36)
        
        push!(yaml_dict["lattice"]["beams"], Dict(
            "charge" => charge_value,
            "brho" => get(brho_values, Int(charge_value), 5.69860938), 
            "matrix" => matrix_flat,
            "centroid" => centroid
        ))
    end
    
    # Save the YAML file
    YAML.write_file(output_filename, yaml_dict)
    
    println("YAML lattice file generated: $output_filename")
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

# Process command line arguments
if length(ARGS) < 1 || length(ARGS) > 2
    println("Usage: julia FLAME_lattice_conversion.jl <input_flame_file> [output_yaml_file]")
    println("       If output file is not specified, it will default to 'jutrack_lattice.yaml'")
    exit(1)
end

input_file = ARGS[1]
output_file = length(ARGS) == 2 ? ARGS[2] : "jutrack_lattice.yaml"

# Check if input file exists
if !isfile(input_file)
    println("Error: Input file '$input_file' does not exist.")
    exit(1)
end

# Run the conversion
println("Converting FLAME lattice file '$input_file' to YAML format...")

try
    success = parse_flame_lattice(input_file, output_file)
    if success
        println("Conversion successful! YAML lattice file saved to '$output_file'")
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