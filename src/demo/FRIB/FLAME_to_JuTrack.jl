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

"""
    parse_flame_lattice(filename::String, output_filename::String)

Parse a FLAME lattice file and convert it to a JuTrack lattice file.
Handles coordinate system transformation from FLAME [mm, rad, mm, rad, rad, MeV/u]
to JuTrack [m, px, m, py, m, dp].
"""
# function parse_flame_lattice(filename::String, output_filename::String)
#     # Read the entire file content
#     file_content = read(filename, String)
    
#     # Initialize collections for elements and global parameters
#     element_defs = Dict{String, Dict{String, Any}}()
#     global_params = Dict{String, Any}()
    
#     # Extract mass number from IonChargeStates
#     mass_number = 0.0  # Will be set based on file content
#     charge_states_pattern = r"IonChargeStates\s*=\s*\[(.*?)\];"
#     charge_match = match(charge_states_pattern, file_content)
    
#     if charge_match !== nothing
#         values = split(charge_match.captures[1], ",")
#         for val_str in values
#             val_str = strip(val_str)
            
#             # Handle division expressions (like 50./124.)
#             if occursin("/", val_str)
#                 parts = split(val_str, "/")
#                 if length(parts) == 2
#                     try
#                         # Extract the denominator which is the mass number
#                         denom_str = replace(strip(parts[2]), r"\.$" => "")
#                         mass_number = parse(Float64, denom_str)
#                         println("Extracted mass number: $mass_number")
#                         break  # We only need to extract it once
#                     catch e
#                         @warn "Could not extract mass number from: $val_str"
#                     end
#                 end
#             end
#         end
#     end
    
#     if mass_number == 0.0
#         # If we couldn't extract from IonChargeStates, look for it in the file
#         mass_pattern = r"mass_number\s*=\s*([0-9\.]+);"
#         mass_match = match(mass_pattern, file_content)
#         if mass_match !== nothing
#             mass_number = parse(Float64, mass_match.captures[1])
#             println("Extracted mass number from explicit declaration: $mass_number")
#         else
#             # Set a reasonable default if we can't find it
#             mass_number = 124.0  # Based on the IonChargeStates from the header
#             println("Using default mass number: $mass_number")
#         end
#     end
    
#     # Extract energy values
#     rest_energy_pattern = r"IonEs\s*=\s*([0-9\.e\+\-]+);"
#     rest_energy_match = match(rest_energy_pattern, file_content)
#     rest_energy_per_nucleon = 931494320.0  # Default value in eV/u
    
#     if rest_energy_match !== nothing
#         try
#             rest_energy_per_nucleon = parse(Float64, rest_energy_match.captures[1])
#             println("Extracted rest energy: $rest_energy_per_nucleon eV/u")
#         catch
#             println("Warning: Could not parse rest energy from file, using default")
#         end
#     end
    
#     kinetic_energy_pattern = r"IonEk\s*=\s*([0-9\.e\+\-]+);"
#     kinetic_energy_match = match(kinetic_energy_pattern, file_content)
#     kinetic_energy_per_nucleon = 227050000.0  # Default value that worked for you
    
#     if kinetic_energy_match !== nothing
#         try
#             kinetic_energy_per_nucleon = parse(Float64, kinetic_energy_match.captures[1])
#             println("Extracted kinetic energy: $kinetic_energy_per_nucleon eV/u")
#         catch
#             println("Warning: Could not parse kinetic energy from file, using default")
#         end
#     end
    
#     # Calculate relativistic parameters
#     gamma = 1.0 + kinetic_energy_per_nucleon / rest_energy_per_nucleon
#     beta = sqrt(1.0 - 1.0 / gamma^2)
#     beta_gamma = beta * gamma
#     println("Relativistic parameters: beta=$beta, gamma=$gamma, beta*gamma=$beta_gamma")
    
#     # Extract RF frequency for phase-to-length conversion
#     rf_freq = 80.5e6  # Default RF frequency in Hz
#     freq_pattern = r"f\s*=\s*([0-9\.e\+\-]+);"
#     for freq_match in eachmatch(freq_pattern, file_content)
#         try
#             rf_freq = parse(Float64, freq_match.captures[1])
#             println("Found RF frequency: $rf_freq Hz")
#             break  # Take the first match
#         catch
#             continue
#         end
#     end
#     rf_wavelength = speed_of_light / rf_freq
#     rf_k = 2 * pi * rf_freq / speed_of_light  # [1/m]
#     println("RF wavelength: $rf_wavelength m")
    
#     # Create the scaling array for coordinate transformation
#     # FLAME: [x(mm), x'(rad), y(mm), y'(rad), φ(rad), dE_k(MeV/u)]
#     # JuTrack: [x(m), px, y(m), py, z(m), dp/p]
#     scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(rest_energy_per_nucleon+kinetic_energy_per_nucleon)]
    
#     println("Coordinate transformation scaling:")
#     display(scaling)
#     println()
    
#     # Extract S matrices for beam envelopes
#     raw_matrices = Dict{String, Matrix{Float64}}()
#     beam_envelopes = Dict{String, Matrix{Float64}}()
#     matrix_names = ["S0", "S1", "S2"]
    
#     for name in matrix_names
#         println("Looking for matrix $name")
#         pattern = Regex("$name\\s*=\\s*\\[(.*?)\\];", "s")
#         m = match(pattern, file_content)
        
#         if m !== nothing
#             println("Found matrix $name")
#             matrix_data = m.captures[1]
            
#             # Create a 6x6 matrix
#             matrix = zeros(Float64, 6, 6)
            
#             # Split into rows
#             rows = split(matrix_data, "\n")
#             rows = [strip(row) for row in rows if !isempty(strip(row))]
            
#             for i in 1:min(6, length(rows))
#                 # Split row into values
#                 values = split(rows[i], ",")
#                 values = [strip(v) for v in values if !isempty(strip(v))]
                
#                 for j in 1:min(6, length(values))
#                     try
#                         # Handle scientific notation and other numeric formats
#                         val_str = values[j]
#                         val = nothing
#                         try
#                             val = parse(Float64, val_str)
#                         catch
#                             try
#                                 val = eval(Meta.parse(val_str))
#                             catch
#                                 println("Warning: Could not parse value '$val_str' at position [$i,$j]")
#                                 val = 0.0
#                             end
#                         end
#                         matrix[i, j] = val
#                     catch e
#                         println("Error parsing value at position [$i,$j]: $e")
#                         matrix[i, j] = 0.0
#                     end
#                 end
#             end
#             lam, u = eigen(matrix)

#             if any(lam .<= 0)
#                     println("Warning: Detected negative eigenvalues in the $name matrix, resetting matrix to zero cross terms")
#                     matrix = zero_non_diagonal_blocks(matrix)
#             end
#             # Store the raw matrix
#             raw_matrices[name] = matrix
            
#             println("Stored raw matrix $name")
#         else
#             println("Warning: Could not find matrix $name")
#         end
#     end
    
#     # Extract and transform BaryCenter vectors
#     barycenters = Dict{String, Vector{Float64}}()
    
#     for i in 0:2
#         name = "BaryCenter$i"
#         pattern = Regex("$name\\s*=\\s*\\[(.*?)\\];", "s")
#         m = match(pattern, file_content)
#         if m !== nothing
#             # Extract the values as a string
#             values_str = m.captures[1]

#             # Split the values string by commas and strip spaces
#             values_str_split = split(values_str, ",")
            
#             # Evaluate each value individually and convert to a number
#             evaluated_values = Float64[]
#             for value in values_str_split
#                 value = strip(value)  # Remove any extra spaces
#                 try
#                     # Evaluate the value if it's a mathematical expression
#                     eval_value = eval(Meta.parse(value))
#                     push!(evaluated_values, eval_value)
#                 catch e
#                     println("Error evaluating value: $value")
#                     push!(evaluated_values, 0.0)  # Use default if parsing fails
#                 end
#             end
            
#             # Remove the 7th element if it exists (the weight/fraction)
#             if length(evaluated_values) >= 7
#                 deleteat!(evaluated_values, 7)
#             end
            
#             # Pad to length 6 if shorter
#             while length(evaluated_values) < 6
#                 push!(evaluated_values, 0.0)
#             end
            
#             # Convert FLAME centroid to JuTrack centroid using the scaling
#             flame_centroid = evaluated_values[1:6]
#             jutrack_centroid = scaling .* flame_centroid
            
#             # Store the transformed centroid
#             barycenters[name] = jutrack_centroid
            
#             println("Transformed $name to JuTrack coordinates: $jutrack_centroid")
#         else
#             println("Warning: Could not find $name")
#         end
#     end
    
#     # Extract charge states and NCharge values
#     charge_states = Float64[]
#     if charge_match !== nothing
#         println("Found IonChargeStates: $(charge_match.captures[1])")
#         values = split(charge_match.captures[1], ",")
#         for val_str in values
#             val_str = strip(val_str)
            
#             # Handle division expressions (like 50./124.)
#             if occursin("/", val_str)
#                 parts = split(val_str, "/")
#                 if length(parts) == 2
#                     try
#                         # Parse numerator as a float
#                         numer = parse(Float64, replace(strip(parts[1]), r"\.$" => ""))
#                         charge_states = push!(charge_states, numer)
#                     catch e
#                         @warn "Could not parse charge state fraction: $val_str"
#                     end
#                 end
#             else
#                 # Try direct parsing
#                 try
#                     push!(charge_states, parse(Float64, val_str))
#                 catch e
#                     @warn "Could not parse charge state: $val_str"
#                 end
#             end
#         end
#         println("Parsed charge states: $charge_states")
#     else
#         println("Warning: Could not find IonChargeStates")
#     end
    
#     # Extract NCharge values
#     ncharge = Float64[]
#     ncharge_pattern = r"NCharge\s*=\s*\[(.*?)\];"
#     ncharge_match = match(ncharge_pattern, file_content)
#     if ncharge_match !== nothing
#         println("Found NCharge: $(ncharge_match.captures[1])")
#         values = split(ncharge_match.captures[1], ",")
#         for val_str in values
#             val_str = strip(val_str)
#             try
#                 push!(ncharge, parse(Float64, val_str))
#             catch e
#                 @warn "Could not parse NCharge value: $val_str"
#             end
#         end
#         println("Parsed NCharge values: $ncharge")
#     else
#         println("Warning: Could not find NCharge")
#     end
    
#     # Extract element definitions and maintain order
#     element_order = String[]
#     element_pattern = r"([^;:]+(?::[^;:]+)?)\s*:\s*([^,;]+)(?:,\s*(.*?))?;"
#     for m in eachmatch(element_pattern, file_content)
#         element_name = String(strip(m[1]))
#         element_type = String(strip(m[2]))
        
#         # Skip LINE definitions, comments, source, or explanatory text
#         if occursin("LINE", element_type) || 
#            element_type == "source" || 
#            startswith(element_name, "#") || 
#            occursin("Note", element_name) ||
#            occursin("Units are", element_name) ||
#            occursin("is:", element_name) ||
#            occursin("are:", element_name) ||
#            length(element_name) > 50  # Likely explanatory text, not an element name
#             continue
#         end
        
#         # Add to the ordered list if not already added
#         if !(element_name in element_order)
#             push!(element_order, element_name)
#         end
        
#         # Parse parameters
#         element_params = Dict{String, Any}()
#         element_params["type"] = element_type
#         element_params["name"] = element_name
        
#         if m[3] !== nothing
#             params_str = m[3]
#             for param in split(params_str, ",")
#                 param = strip(param)
#                 if occursin("=", param)
#                     param_parts = split(param, "=", limit=2)
#                     param_name = String(strip(param_parts[1]))
#                     param_value_str = String(strip(param_parts[2]))
                    
#                     if param_name == "bg"  # Skip beam gamma parameter
#                         continue
#                     end

#                     # Try to convert to appropriate type
#                     param_value = try
#                         parse(Float64, param_value_str)
#                     catch
#                         # Remove quotes and try numeric evaluation
#                         clean_value = replace(param_value_str, "\"" => "")
#                         try
#                             if occursin(r"[0-9+\-*/\^\.e]", clean_value) && !occursin(r"[a-zA-Z]", clean_value)
#                                 eval(Meta.parse(clean_value))
#                             else
#                                 clean_value
#                             end
#                         catch
#                             clean_value
#                         end
#                     end
                    
#                     element_params[param_name] = param_value
#                 end
#             end
#         end
        
#         element_defs[element_name] = element_params
#     end
    
#     # Extract energy and mass values with proper unit conversion
#     total_kinetic_energy = kinetic_energy_per_nucleon * mass_number  # Total kinetic energy in eV
#     total_rest_energy = rest_energy_per_nucleon * mass_number  # Total rest energy in eV
#     total_energy = total_kinetic_energy + total_rest_energy  # Total energy in eV
    
#     # Now generate the JuTrack lattice file
#     open(output_filename, "w") do file
#         # Write header
#         write(file, """
#         # JuTrack lattice file converted from FLAME format
#         # Original file: $filename
#         # Generated on: $(Dates.now())
#         # Mass number: $mass_number
#         # Coordinate system conversion: FLAME [mm, rad, mm, rad, rad, MeV/u] -> JuTrack [m, px, m, py, m, dp]
#         # Relativistic factors: beta=$beta, gamma=$gamma, beta*gamma=$beta_gamma
#         # RF frequency: $rf_freq Hz, wavelength: $rf_wavelength m
        
#         include("../../JuTrack.jl")
#         using .JuTrack
#         using LinearAlgebra
#         using Random
#         using Distributions
#         using CairoMakie
#         using Statistics
#         using ProgressMeter
        
#         """
#         )
        
#         # Write function to create the lattice
#         write(file, """
#         function create_lattice()
#             lattice = []
        
#             # Define elements
#         """)
        
#         # Define all elements in order
#         for element_name in element_order
#             if haskey(element_defs, element_name)
#                 params = element_defs[element_name]
                
#                 # Generate a safe name for Julia
#                 safe_name = replace(element_name, r"[: #\.]" => "_")
#                 safe_name = replace(safe_name, r"[\(\)]" => "")
                
#                 # Map FLAME element type to JuTrack type
#                 jutrack_type = get_jutrack_type(params["type"])
                
#                 # Generate parameter list
#                 param_str = "name=\"$(element_name)\""
                
#                 for (param_name, param_value) in params
#                     if param_name != "type" && param_name != "name"
#                         # Map parameter names and apply unit conversions if needed
#                         jutrack_param = get_jutrack_param(param_name, param_value, params["type"], mass_number)
#                         if !isnothing(jutrack_param)
#                             param_str *= ", $(jutrack_param)"
#                         end
#                     end
#                 end
                
#                 # Create element
#                 write(file, "    $safe_name = $jutrack_type($param_str)\n")
                
#                 # Add element to the lattice
#                 write(file, "    push!(lattice, $safe_name)\n")
#             end
#         end
        
#         # Write raw matrices
#         write(file, "\n    # Define FLAME raw matrices\n")
#         for (matrix_name, matrix) in raw_matrices
#             write(file, "    $(matrix_name)_raw = [\n")
#             for i in 1:6
#                 write(file, "        ")
#                 for j in 1:6
#                     write(file, "$(matrix[i, j]) ")
#                 end
#                 write(file, ";\n")
#             end
#             write(file, "    ]\n\n")
#         end
        
#         # Write the RF and beam parameters
#         write(file, """
#             # Define RF and beam parameters
#             rf_freq = $rf_freq  # RF frequency in Hz
#             rf_k = 2 * pi * rf_freq / 299792458.0  # [1/m]
#             IonEs = $rest_energy_per_nucleon  # Nucleon mass [eV/u]
#             IonEk = $kinetic_energy_per_nucleon  # Kinetic energy [eV/u]
#             gamma = 1.0 + IonEk / IonEs  # Lorentz factor
#             beta = sqrt(1.0 - 1.0 / (gamma * gamma))  # Velocity factor
            
#             # Define scaling factors for coordinate transformation
#             scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(IonEs+IonEk)]
            
#             # Transform FLAME matrices to JuTrack matrices
#         """)
        
#         # Write code to transform matrices
#         for i in 0:2
#             matrix_name = "S$i"
#             write(file, """
#                 # Transform $matrix_name
#                 $(matrix_name)_matrix = diagm(scaling) * (($(matrix_name)_raw + $(matrix_name)_raw') ./ 2) * diagm(scaling)'  # Apply scaling and symmetrize
            
#             """)
#         end
        
#         # Write beam creation code with exact approach
#         write(file, """
#             # Create beam objects for each charge state
#         """)
        
#         for i in 1:min(3, length(charge_states))
#             # Get barycenter data
#             bary_key = "BaryCenter$(i-1)"
#             centroid = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Default centroid
            
#             if haskey(barycenters, bary_key)
#                 centroid = barycenters[bary_key]
#             end
            
#             # Get charge value and weight
#             charge_value = charge_states[i]
#             weight_value = (i <= length(ncharge)) ? ncharge[i] : 1.0
            
#             write(file, """
#                 # Generate beam$i using eigendecomposition
#                 Random.seed!($(42+i-1))  # Ensure reproducibility
                
#                 # Get eigendecomposition of the transformed matrix
#                 lam$i, u$i = eigen(S$(i-1)_matrix)

#                 # Generate random particles 
#                 nparticles = 10000000
#                 dis$i = Matrix{Float64}(undef, nparticles, 6)
#                 for d in 1:6
#                     dis$i[:, d] .= randn(nparticles)
#                 end
                
#                 # Calculate the initial 2nd moment matrix 
#                 moment2nd$i = zeros(6, 6)
#                 for d1 in 1:6
#                     for d2 in 1:6
#                         moment2nd$i[d1, d2] = mean(dis$i[:, d1] .* dis$i[:, d2])
#                     end
#                 end

#                 lam$i,u$i = eigen(moment2nd$i)
#                 transformation$i = u$i*diagm(1.0 ./ sqrt.(lam$i)) *u$i'
#                 dis$i = dis$i * transformation$i'
                
#                 lam$i,u$i = eigen(S$(i-1)_matrix)
#                 # Create transformation matrix from eigendecomposition
#                 transformation$i = u$i * Diagonal(sqrt.(lam$i)) * u$i'
                
#                 # Apply transformation to get desired covariance
#                 dis$i = dis$i * transformation$i'
                
#                 # Calculate the final 2nd moment matrix to verify transformation
#                 for d1 in 1:6
#                     for d2 in 1:6
#                         moment2nd$i[d1, d2] = mean(dis$i[:, d1] .* dis$i[:, d2])
#                     end
#                 end
                
#                 # Calculate and display the ratio to verify accuracy
#                 # ratio$i = moment2nd$i ./ S$(i-1)_matrix
#                 ratio$i = ifelse.((moment2nd$i .<= 1e-20) .& (S$(i-1)_matrix .== 0.0), 1.0, moment2nd$i ./ S$(i-1)_matrix)
#                 println("Beam $i second moment ratio:")
#                 display(ratio$i)
                
#                 # Add centroid
#                 dis$i .+= reshape([$(centroid[1]), $(centroid[2]), $(centroid[3]), $(centroid[4]), $(centroid[5]), $(centroid[6])], 1, 6)
                
#                 # Create beam with the generated particles - using a subset for efficiency
#                 beam$i = Beam(r=dis$i[1:10000,:], np=10000, energy=$total_kinetic_energy, charge=$charge_value, mass=$total_rest_energy)
#                 get_centroid!(beam$i)
#                 get_emittance!(beam$i)
            
#             """)
#         end
        
#         write(file, "    return lattice, beam1, beam2, beam3\n")
#         write(file, "end\n\n")
        
#         # Rest of the file remains the same
#         write(file, """
#         # Run the simulation

#         function propagate_beam(lattice, beam, np)
#             # Create expanded lattice
#             expanded_lattice = []
#             expanded_indices = []  # Track which original element each expanded element belongs to
            
#             @showprogress for (i, element) in enumerate(lattice)
#                 if element.len > 0
#                     # Calculate segments needed
#                     segments_by_count = 10
#                     segments_by_length = ceil(Int, element.len / 0.1)
#                     num_segments = max(segments_by_count, segments_by_length)
#                     segment_length = element.len / num_segments
                    
#                     # Create smaller elements
#                     for j in 1:num_segments
#                         small_element = deepcopy(element)
#                         small_element.len = segment_length
#                         push!(expanded_lattice, small_element)
#                         push!(expanded_indices, i)  # Record original index
#                     end
#                 else
#                     # Zero-length elements
#                     push!(expanded_lattice, element)
#                     push!(expanded_indices, i)
#                 end
#             end
            
#             # Propagate through expanded lattice
#             n_expanded = length(expanded_lattice)
#             expanded_beam_rms = zeros(n_expanded, 6)
#             expanded_floor_distance = zeros(n_expanded)
#             expanded_twi = zeros(n_expanded, 9)
#             flat_particles = zeros(6 * np)
            
#             # Propagate through expanded lattice
#             for i in eachindex(expanded_lattice)
#                 expanded_twi[i, :] .= twiss_beam(beam)
#                 flat_particles .= collect(Iterators.flatten(eachrow(beam.r)))
                
#                 # Pass particles through element
#                 pass!(expanded_lattice[i], flat_particles, np, beam)
#                 beam.r = reshape(flat_particles, 6, np)'
                
#                 # Store RMS values
#                 for dim in 1:6
#                     expanded_beam_rms[i, dim] = sqrt(mean(beam.r[:,dim].^2))
#                 end
                
#                 # Calculate floor distance
#                 if i == 1
#                     expanded_floor_distance[i] = expanded_lattice[i].len
#                 else
#                     expanded_floor_distance[i] = expanded_floor_distance[i-1] + expanded_lattice[i].len
#                 end
#             end
            
#             # Return everything needed
#             return expanded_beam_rms, expanded_floor_distance, beam, expanded_twi, 
#                 lattice, expanded_lattice, expanded_indices
#         end

#         function run_simulation()
#             lattice, beam1, beam2, beam3 = create_lattice()
            
#             beam1_rms, floor_distance, beam1, twi1, original_lattice, _, expanded_indices1 = 
#                 propagate_beam(lattice, beam1, beam1.np)
            
#             beam2_rms, _, beam2, twi2, _, _, expanded_indices2 = 
#                 propagate_beam(lattice, beam2, beam2.np)
            
#             beam3_rms, _, beam3, twi3, _, _, expanded_indices3 = 
#                 propagate_beam(lattice, beam3, beam3.np)

#             return original_lattice, beam1, beam2, beam3, 
#                 beam1_rms, beam2_rms, beam3_rms, 
#                 twi1, twi2, twi3, floor_distance, expanded_indices1
#         end

#         # plotting wrapper function
#         function plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, 
#                                             twi1, twi2, twi3, floor_distance,
#                                             original_lattice, expanded_indices, 
#                                             combined=false)
#             # Filter expanded data up to the specified original element
#             mask = expanded_indices .<= end_ele
            
#             # Use the filtered data for plots but original lattice for floor layout
#             return plot_multibeam_data(
#                 floor_distance[mask], 
#                 beam1_rms[mask,:], 
#                 beam2_rms[mask,:], 
#                 beam3_rms[mask,:], 
#                 twi1[mask,:], 
#                 twi2[mask,:], 
#                 twi3[mask,:], 
#                 original_lattice[1:end_ele],  # Use original lattice for floor plot
#                 combined=combined
#             )
#         end

#         function visualize_beam_properties(beam::Beam)
#             # Create a multi-panel figure
#             fig = Figure(size=(900, 600))
            
#             # Phase space plots (x-px, y-py, z-dp)
#             ax1 = Axis(fig[1, 1], title="x-px phase space", xlabel="x [mm]", ylabel="px")
#             # Convert m to mm for x-axis
#             scatter!(ax1, beam.r[:,1] .* 1000, beam.r[:,2], markersize=1, strokewidth=0, alpha=0.5)
            
#             ax2 = Axis(fig[1, 2], title="y-py phase space", xlabel="y [mm]", ylabel="py")
#             # Convert m to mm for y-axis
#             scatter!(ax2, beam.r[:,3] .* 1000, beam.r[:,4], markersize=1, strokewidth=0, alpha=0.5)
            
#             ax3 = Axis(fig[1, 3], title="z-dp phase space", xlabel="z [mm]", ylabel="dp/p")
#             # Convert m to mm for z-axis
#             scatter!(ax3, beam.r[:,5] .* 1000, beam.r[:,6], markersize=1, strokewidth=0, alpha=0.5)
            
#             # Projections (histograms)
#             ax4 = Axis(fig[2, 1], title="x distribution", xlabel="x [mm]")
#             # Convert m to mm for x-axis
#             hist!(ax4, beam.r[:,1] .* 1000, bins=50, color=(:blue, 0.7))
            
#             ax5 = Axis(fig[2, 2], title="y distribution", xlabel="y [mm]")
#             # Convert m to mm for y-axis
#             hist!(ax5, beam.r[:,3] .* 1000, bins=50, color=(:blue, 0.7))
            
#             ax6 = Axis(fig[2, 3], title="z distribution", xlabel="z [mm]")
#             # Convert m to mm for z-axis
#             hist!(ax6, beam.r[:,5] .* 1000, bins=50, color=(:blue, 0.7))
            
#             # Return the figure
#             return fig
#         end

#         function plot_multibeam_data(floor_length, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, lattice; combined = true)
#             if combined != true    
#                 # Define beam colors and labels (shared across all plots)
#                 beam_colors = [:blue, :orange, :green]
#                 beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]
                
#                 # Define element colors and heights (shared across all plots)
#                 element_colors = Dict(
#                     "SBEND" => :purple,
#                     "KQUAD" => :red,
#                     "ORBTRIM" => :green,
#                     "MARKER" => :gray,
#                     "DRIFT" => :white,
#                     "default" => :gray
#                 )
                
#                 element_heights = Dict(
#                     "SBEND" => 0.8,
#                     "KQUAD" => 0.7,
#                     "ORBTRIM" => 1.0,
#                     "MARKER" => 0.4,
#                     "DRIFT" => 0.3,
#                     "default" => 0.5
#                 )
                
#                 # Helper function to create a standard plot structure
#                 function create_plot(top_data1, top_data2, top_data3, bottom_data1, bottom_data2, bottom_data3, 
#                                     top_label, bottom_label, title)
                    
#                     f = Figure(size = (1200, 400))
#                     gl = f[1, 1] = GridLayout()
#                     legend_layout = f[1, 2] = GridLayout()
                    
#                     # Top plot
#                     ax_top = Axis(gl[1, 1], xlabel = "", ylabel = top_label)
#                     lines!(ax_top, floor_length, top_data1, color = beam_colors[1], linewidth = 2)
#                     lines!(ax_top, floor_length, top_data2, color = beam_colors[2], linewidth = 2)
#                     lines!(ax_top, floor_length, top_data3, color = beam_colors[3], linewidth = 2)
#                     hidexdecorations!(ax_top)
                    
#                     # Create floorline plot in the middle
#                     floor_ax = Axis(gl[2, 1], 
#                                 xlabel = "",
#                                 ylabel = "",
#                                 yticklabelsvisible = false,
#                                 yticksvisible = false)
                    
#                     # Track element types for legend
#                     element_types = Dict{String, Any}()
                    
#                     # Draw the floorline
#                     curr_pos = 0.0
#                     for ele in lattice
#                         # Get element type
#                         ele_type = string(typeof(ele).name.name)
                        
#                         # Determine color and height
#                         color = get(element_colors, ele_type, element_colors["default"])
#                         height = get(element_heights, ele_type, element_heights["default"])
                        
#                         # Make sure elements don't overlap by using a minimum width
#                         ele_width = ele.len #max(ele.len, 0.05)
                        
#                         # Draw rectangle for element based on element type and properties
#                         if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
#                             # Use the original quad height
#                             quad_height = element_heights["KQUAD"]
#                             drift_height = element_heights["DRIFT"]
                            
#                             # Calculate offset to align with other elements
#                             # This will make it overlap the x-axis slightly to maintain alignment
#                             if ele.k1 > 0
#                                 # For focusing quads: align bottom with the bottom of drift elements
#                                 y_pos = -drift_height/2
#                             else
#                                 # For defocusing quads: align top with the top of drift elements
#                                 y_pos = drift_height/2 - quad_height
#                             end
                            
#                             rect = Rect(curr_pos, y_pos, ele_width, quad_height)
#                         elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
#                             # Use the original sextupole height
#                             sext_height = get(element_heights, "KSEXT", element_heights["default"])
#                             drift_height = element_heights["DRIFT"]
                            
#                             # Calculate offset to align with other elements
#                             if ele.k2 > 0
#                                 # For focusing sextupoles: align bottom with the bottom of drift elements
#                                 y_pos = -drift_height/2
#                             else
#                                 # For defocusing sextupoles: align top with the top of drift elements
#                                 y_pos = drift_height/2 - sext_height
#                             end
                            
#                             rect = Rect(curr_pos, y_pos, ele_width, sext_height)
#                         else
#                             # All other elements centered on x-axis as before
#                             rect = Rect(curr_pos, -height/2, ele_width, height)
#                         end
                        
#                         element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                        
#                         # Track for legend (only if not already tracked)
#                         if !haskey(element_types, ele_type)
#                             element_types[ele_type] = element
#                         end
                        
#                         # Update position - ensure we advance by at least the element width
#                         # to avoid overlapping due to position calculations
#                         curr_pos += ele_width
#                     end
                    
#                     hidexdecorations!(floor_ax)
                    
#                     # Bottom plot
#                     ax_bottom = Axis(gl[3, 1], xlabel = "z [m]", ylabel = bottom_label, yreversed = true)
#                     lines!(ax_bottom, floor_length, bottom_data1, color = beam_colors[1], linewidth = 2)
#                     lines!(ax_bottom, floor_length, bottom_data2, color = beam_colors[2], linewidth = 2)
#                     lines!(ax_bottom, floor_length, bottom_data3, color = beam_colors[3], linewidth = 2)
                    
#                     # Link x axes
#                     linkxaxes!(ax_top, floor_ax, ax_bottom)
                    
#                     # Add beam legend
#                     legend_entries = [
#                         LineElement(color = c, linewidth = 2) for c in beam_colors
#                     ]
                    
#                     leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                    
#                     # Add element legend (filtering out MARKER and DRIFT)
#                     element_entries = []
#                     element_names = []
#                     sorted_names = sort(collect(keys(element_types)))
#                     for name in sorted_names
#                         if !(name in ["MARKER", "DRIFT"])
#                             element = element_types[name]
#                             push!(element_entries, element)
#                             push!(element_names, name)
#                         end
#                     end
                    
#                     if !isempty(element_entries)
#                         element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
#                     end
                    
#                     # Adjust spacing between legends
#                     rowgap!(legend_layout, 5)
                    
#                     # Set row sizes
#                     rowsize!(gl, 1, 100)  # Top plot
#                     rowsize!(gl, 2, 30)   # Elements plot (smaller)
#                     rowsize!(gl, 3, 100)  # Bottom plot
                    
#                     # Set title
#                     # f.title = title
                    
#                     return f
#                 end
                
#                 # Create the four plots
#                 rms_plot = create_plot(
#                     beam1_rms[:,1] .* 1e3,
#                     beam2_rms[:,1] .* 1e3,
#                     beam3_rms[:,1] .* 1e3,
#                     beam1_rms[:,3] .* 1e3,
#                     beam2_rms[:,3] .* 1e3,
#                     beam3_rms[:,3] .* 1e3,
#                     "RMS X [mm]",
#                     "RMS Y [mm]",
#                     "RMS Plot"
#                 )
                
#                 beta_plot = create_plot(
#                     twi1[:,1],
#                     twi2[:,1],
#                     twi3[:,1],
#                     twi1[:,4],
#                     twi2[:,4],
#                     twi3[:,4],
#                     L"\\beta_x",
#                     L"\\beta_y",
#                     "Beta Functions"
#                 )
                
#                 alpha_plot = create_plot(
#                     twi1[:,2],
#                     twi2[:,2],
#                     twi3[:,2],
#                     twi1[:,5],
#                     twi2[:,5],
#                     twi3[:,5],
#                     L"\\alpha_x",
#                     L"\\alpha_y",
#                     "Alpha Functions"
#                 )
                
#                 emittance_plot = create_plot(
#                     twi1[:,3],
#                     twi2[:,3],
#                     twi3[:,3],
#                     twi1[:,6],
#                     twi2[:,6],
#                     twi3[:,6],
#                     L"\\epsilon_x",
#                     L"\\epsilon_y",
#                     "Emittance"
#                 )
                
#                 # Return all plots
#                 return Dict(
#                     "rms" => rms_plot,
#                     "beta" => beta_plot,
#                     "alpha" => alpha_plot,
#                     "emittance" => emittance_plot
#                 )
#             else
#                 # Create one large figure
#                 f = Figure(size = (1200, 1200))
                
#                 # Create the main layout
#                 gl = f[1, 1] = GridLayout()
#                 legend_layout = f[1, 2] = GridLayout()
                
#                 # Define beam colors and labels
#                 beam_colors = [:blue, :orange, :green]
#                 beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]
                
#                 # Define element colors and heights
#                 element_colors = Dict(
#                     "SBEND" => :purple,
#                     "KQUAD" => :red,
#                     "ORBTRIM" => :green,
#                     "MARKER" => :gray,
#                     "DRIFT" => :white,
#                     "default" => :gray
#                 )
                
#                 element_heights = Dict(
#                     "SBEND" => 0.8,
#                     "KQUAD" => 0.7,
#                     "ORBTRIM" => 0.6,
#                     "MARKER" => 0.4,
#                     "DRIFT" => 0.3,
#                     "default" => 0.5
#                 )
                
#                 # Create plots in sequence (8 total + 1 lattice in the middle)
#                 # 1. RMS X
#                 ax_rms_x = Axis(gl[1, 1], xlabel = "", ylabel = "RMS X [mm]", title = "RMS Values")
#                 lines!(ax_rms_x, floor_length, beam1_rms[:,1] .* 1e3, color = beam_colors[1], linewidth = 2)
#                 lines!(ax_rms_x, floor_length, beam2_rms[:,1] .* 1e3, color = beam_colors[2], linewidth = 2)
#                 lines!(ax_rms_x, floor_length, beam3_rms[:,1] .* 1e3, color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_rms_x)
                
#                 # 2. RMS Y
#                 ax_rms_y = Axis(gl[2, 1], xlabel = "", ylabel = "RMS Y [mm]")
#                 lines!(ax_rms_y, floor_length, beam1_rms[:,3] .* 1e3, color = beam_colors[1], linewidth = 2)
#                 lines!(ax_rms_y, floor_length, beam2_rms[:,3] .* 1e3, color = beam_colors[2], linewidth = 2)
#                 lines!(ax_rms_y, floor_length, beam3_rms[:,3] .* 1e3, color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_rms_y)
                
#                 # 3. Beta X
#                 ax_beta_x = Axis(gl[3, 1], xlabel = "", ylabel = L"\\beta_x", title = "Beta Functions")
#                 lines!(ax_beta_x, floor_length, twi1[:,1], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_beta_x, floor_length, twi2[:,1], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_beta_x, floor_length, twi3[:,1], color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_beta_x)
                
#                 # 4. Beta Y
#                 ax_beta_y = Axis(gl[4, 1], xlabel = "", ylabel = L"\\beta_y")
#                 lines!(ax_beta_y, floor_length, twi1[:,4], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_beta_y, floor_length, twi2[:,4], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_beta_y, floor_length, twi3[:,4], color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_beta_y)
                
#                 # 5. Lattice plot in the middle
#                 floor_ax = Axis(gl[5, 1], 
#                             xlabel = "",
#                             ylabel = "",
#                             yticklabelsvisible = false,
#                             yticksvisible = false)
                
#                 # Track element types for legend
#                 element_types = Dict{String, Any}()
                
#                 # Draw the floorline
#                 curr_pos = 0.0
#                 for ele in lattice
#                     # Get element type
#                     ele_type = string(typeof(ele).name.name)
                    
#                     # Determine color and height
#                     color = get(element_colors, ele_type, element_colors["default"])
#                     height = get(element_heights, ele_type, element_heights["default"])
                    
#                     # Use a minimum width to prevent tiny elements and avoid overlap
#                     ele_width = ele.len #max(ele.len, 0.05)
                    
#                     # Draw rectangle for element based on element type and properties
#                     if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
#                         # Use the original quad height, not drift height
#                         quad_height = element_heights["KQUAD"]
#                         if ele.k1 > 0
#                             # Position touching x-axis and extending upward with full quad height
#                             rect = Rect(curr_pos, 0, ele_width, quad_height)
#                         else
#                             # Position touching x-axis and extending downward with full quad height
#                             rect = Rect(curr_pos, -quad_height, ele_width, quad_height)
#                         end
#                     elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
#                         # Use the original sextupole height
#                         sext_height = get(element_heights, "KSEXT", element_heights["default"])
#                         if ele.k2 > 0
#                             # Position touching x-axis and extending upward with full sextupole height
#                             rect = Rect(curr_pos, 0, ele_width, sext_height)
#                         else
#                             # Position touching x-axis and extending downward with full sextupole height
#                             rect = Rect(curr_pos, -sext_height, ele_width, sext_height)
#                         end
#                     else
#                         # All other elements centered on x-axis as before
#                         rect = Rect(curr_pos, -height/2, ele_width, height)
#                     end
                    
#                     element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                    
#                     # Track for legend (only if not already tracked)
#                     if !haskey(element_types, ele_type)
#                         element_types[ele_type] = element
#                     end
                    
#                     # Update position - use the element width to ensure no overlap
#                     curr_pos += ele_width
#                 end
#                 hidexdecorations!(floor_ax)
                
#                 # 6. Alpha X
#                 ax_alpha_x = Axis(gl[6, 1], xlabel = "", ylabel = L"\\alpha_x", title = "Alpha Functions")
#                 lines!(ax_alpha_x, floor_length, twi1[:,2], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_alpha_x, floor_length, twi2[:,2], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_alpha_x, floor_length, twi3[:,2], color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_alpha_x)
                
#                 # 7. Alpha Y
#                 ax_alpha_y = Axis(gl[7, 1], xlabel = "", ylabel = L"\\alpha_y")
#                 lines!(ax_alpha_y, floor_length, twi1[:,5], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_alpha_y, floor_length, twi2[:,5], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_alpha_y, floor_length, twi3[:,5], color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_alpha_y)
                
#                 # 8. Emittance X
#                 ax_emit_x = Axis(gl[8, 1], xlabel = "", ylabel = L"\\epsilon_x", title = "Emittance")
#                 lines!(ax_emit_x, floor_length, twi1[:,3], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_emit_x, floor_length, twi2[:,3], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_emit_x, floor_length, twi3[:,3], color = beam_colors[3], linewidth = 2)
#                 hidexdecorations!(ax_emit_x)
                
#                 # 9. Emittance Y (bottom plot)
#                 ax_emit_y = Axis(gl[9, 1], xlabel = "z [m]", ylabel = L"\\epsilon_y")
#                 lines!(ax_emit_y, floor_length, twi1[:,6], color = beam_colors[1], linewidth = 2)
#                 lines!(ax_emit_y, floor_length, twi2[:,6], color = beam_colors[2], linewidth = 2)
#                 lines!(ax_emit_y, floor_length, twi3[:,6], color = beam_colors[3], linewidth = 2)
                
#                 # Link all x axes
#                 for ax in [ax_rms_x, ax_rms_y, ax_beta_x, ax_beta_y, floor_ax, 
#                         ax_alpha_x, ax_alpha_y, ax_emit_x, ax_emit_y]
#                     linkxaxes!(ax, ax_rms_x)
#                 end
                
#                 # Add beam legend
#                 legend_entries = [
#                     LineElement(color = c, linewidth = 2) for c in beam_colors
#                 ]
                
#                 leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                
#                 # Add element legend (excluding MARKER and DRIFT)
#                 element_entries = []
#                 element_names = []
#                 sorted_names = sort(collect(keys(element_types)))
#                 for name in sorted_names
#                     if !(name in ["MARKER", "DRIFT"])
#                         element = element_types[name]
#                         push!(element_entries, element)
#                         push!(element_names, name)
#                     end
#                 end
                
#                 if !isempty(element_entries)
#                     element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
#                 end
                
#                 # Adjust spacing between legends
#                 rowgap!(legend_layout, 5)
                
#                 # Set row sizes (make lattice row smaller)
#                 rowsize!(gl, 5, 30)  # Lattice plot is smaller
                
                
#                 return f
#             end
#         end

#         function reconstruct_original_floor_distance(expanded, distance)
#             # Find the maximum original index
#             max_index = maximum(expanded)
            
#             # Initialize the original vector with zeros
#             original = zeros(eltype(distance), max_index)
            

#             for i in 1:length(expanded)
#                 original[expanded[i]] = distance[i]
#             end
            
#             return original
#         end

#         function FLAME_to_Julia_matrix(matrix)
#             new_mat = reshape(matrix, 7, 7)
#             return new_mat[1:end-1, 1:end-1]
#         end

#         function zero_non_diagonal_blocks(matrix)
#             # Check that we have a square matrix
#             rows, cols = size(matrix)
#             if rows != cols || rows % 2 != 0
#                 error("Matrix must be square with even dimensions")
#             end
            
#             # Create a copy of the matrix
#             result = copy(matrix)
            
#             # Size of sub-matrices
#             block_size = 2
            
#             # Zero out all elements not in diagonal blocks
#             for i in 1:rows
#                 for j in 1:cols
#                     # Calculate which block this element belongs to
#                     block_i = ceil(Int, i/block_size)
#                     block_j = ceil(Int, j/block_size)
                    
#                     # If not in a diagonal block, set to zero
#                     if block_i != block_j
#                         result[i, j] = 0
#                     end
#                 end
#             end
            
#             return result
#         end
        
#         # Create the lattice and beams
#         lattice, beam1, beam2, beam3 = create_lattice();
        

#         begin
#             p1 = visualize_beam_properties(beam1);
#             p2 = visualize_beam_properties(beam2);
#             p3 = visualize_beam_properties(beam3);
#             display(p1)
#             display(p2)
#             display(p3)
#         end;


#         # Run simulation
#         lattice, beam1, beam2, beam3, beam1_rms, beam2_rms, beam3_rms, 
#         twi1, twi2, twi3, floor_distance, expanded_indices = run_simulation();

#         # Plot using original indexing
#         end_ele = 200;
#         plots_jutrack = plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, false);
#         plots_jutrack["rms"]
#         plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, true)

#         """)
#     end
    
#     println("JuTrack lattice file generated: $output_filename")
#     return true
# end

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
    kinetic_energy_per_nucleon = 227050000.0  # Default value that worked for you
    
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
            
            # Handle division expressions (like 50./124.)
            if occursin("/", val_str)
                parts = split(val_str, "/")
                if length(parts) == 2
                    try
                        # Extract numerator as a float
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
                    
                    # Instead of skipping bg, store it in element_params
                    # if param_name == "bg"  # Skip beam gamma parameter
                    #     continue
                    # end

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
        using Statistics
        using ProgressMeter
        
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
                
                # Get bg value if it exists for this element
                local element_bg = get(params, "bg", nothing)
                
                for (param_name, param_value) in params
                    if param_name != "type" && param_name != "name" && param_name != "bg"
                        # Map parameter names and apply unit conversions if needed
                        jutrack_param = get_jutrack_param(param_name, param_value, params["type"], mass_number)
                        
                        # if params["type"] == "orbtrim" && param_name == "tm_xkick"
                        #     scaled_value = param_value / 5.69860938 # Divided by its constant Bρ value
                        #     jutrack_param = get_jutrack_param(param_name, scaled_value, params["type"], mass_number)
                        # end
                        # if params["type"] == "orbtrim" && param_name == "tm_ykick"
                        #     scaled_value = param_value / 5.69860938 # Divided by its constant Bρ value
                        #     jutrack_param = get_jutrack_param(param_name, scaled_value, params["type"], mass_number)
                        # end

                        if params["type"] == "quadrupole" && param_name == "B2"
                            scaled_value = param_value / 5.69860938 # Divided by its constant Bρ value
                            jutrack_param = get_jutrack_param(param_name, scaled_value, params["type"], mass_number)
                            println("Scaling k1 for $(element_name): $param_value -> $scaled_value ")
                        end
                        # For SBEND elements, scale the bending angle (phi) by the ratio of actual beta_gamma to reference bg
                        if params["type"] == "sbend" && element_bg !== nothing && 
                           (param_name == "phi" || param_name == "angle")
                            scaled_value = param_value * π / 180.0
                            jutrack_param = get_jutrack_param(param_name, scaled_value, params["type"], mass_number)
                            println("Scaling bending angle for $(element_name): $param_value -> $scaled_value (reference_bg=$element_bg, actual_bg=$beta_gamma)")
                        end
                        
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
        
        # Write raw matrices
        write(file, "\n    # Define FLAME raw matrices\n")
        for (matrix_name, matrix) in raw_matrices
            write(file, "    $(matrix_name)_raw = [\n")
            for i in 1:6
                write(file, "        ")
                for j in 1:6
                    write(file, "$(matrix[i, j]) ")
                end
                write(file, ";\n")
            end
            write(file, "    ]\n\n")
        end
        
        # Write the RF and beam parameters
        write(file, """
            # Define RF and beam parameters
            rf_freq = $rf_freq  # RF frequency in Hz
            rf_k = 2 * pi * rf_freq / 299792458.0  # [1/m]
            IonEs = $rest_energy_per_nucleon  # Nucleon mass [eV/u]
            IonEk = $kinetic_energy_per_nucleon  # Kinetic energy [eV/u]
            gamma = 1.0 + IonEk / IonEs  # Lorentz factor
            beta = sqrt(1.0 - 1.0 / (gamma * gamma))  # Velocity factor
            
            # Define scaling factors for coordinate transformation
            scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(IonEs+IonEk)]
            
            # Transform FLAME matrices to JuTrack matrices
        """)
        
        # Write code to transform matrices
        for i in 0:2
            matrix_name = "S$i"
            write(file, """
                # Transform $matrix_name
                $(matrix_name)_matrix = diagm(scaling) * (($(matrix_name)_raw + $(matrix_name)_raw') ./ 2) * diagm(scaling)'  # Apply scaling and symmetrize
            
            """)
        end
        
        # Write beam creation code with exact approach
        write(file, """
            # Create beam objects for each charge state
        """)
        
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
            
            write(file, """
                # Generate beam$i using eigendecomposition
                Random.seed!($(42+i-1))  # Ensure reproducibility
                
                # Get eigendecomposition of the transformed matrix
                lam$i, u$i = eigen(S$(i-1)_matrix)

                # Generate random particles 
                nparticles = 10000000
                dis$i = Matrix{Float64}(undef, nparticles, 6)
                for d in 1:6
                    dis$i[:, d] .= randn(nparticles)
                end
                
                # Calculate the initial 2nd moment matrix 
                moment2nd$i = zeros(6, 6)
                for d1 in 1:6
                    for d2 in 1:6
                        moment2nd$i[d1, d2] = mean(dis$i[:, d1] .* dis$i[:, d2])
                    end
                end

                lam$i,u$i = eigen(moment2nd$i)
                transformation$i = u$i*diagm(1.0 ./ sqrt.(lam$i)) *u$i'
                dis$i = dis$i * transformation$i'
                
                lam$i,u$i = eigen(S$(i-1)_matrix)
                # Create transformation matrix from eigendecomposition
                transformation$i = u$i * Diagonal(sqrt.(lam$i)) * u$i'
                
                # Apply transformation to get desired covariance
                dis$i = dis$i * transformation$i'
                
                # Calculate the final 2nd moment matrix to verify transformation
                for d1 in 1:6
                    for d2 in 1:6
                        moment2nd$i[d1, d2] = mean(dis$i[:, d1] .* dis$i[:, d2])
                    end
                end
                
                # Calculate and display the ratio to verify accuracy
                # ratio$i = moment2nd$i ./ S$(i-1)_matrix
                ratio$i = ifelse.((moment2nd$i .<= 1e-20) .& (S$(i-1)_matrix .== 0.0), 1.0, moment2nd$i ./ S$(i-1)_matrix)
                println("Beam $i second moment ratio:")
                display(ratio$i)
                
                # Add centroid
                dis$i .+= reshape([$(centroid[1]), $(centroid[2]), $(centroid[3]), $(centroid[4]), $(centroid[5]), $(centroid[6])], 1, 6)
                
                # Create beam with the generated particles - using a subset for efficiency
                beam$i = Beam(r=dis$i[1:10000,:], np=10000, energy=$total_kinetic_energy, charge=$charge_value, mass=$total_rest_energy)
                get_centroid!(beam$i)
                get_emittance!(beam$i)
            
            """)
        end
        
        write(file, "    return lattice, beam1, beam2, beam3\n")
        write(file, "end\n\n")
        
        # Rest of the file remains the same
        write(file, """
        # Run the simulation

        function propagate_beam(lattice, beam, np)
            # Create expanded lattice
            expanded_lattice = []
            expanded_indices = []  # Track which original element each expanded element belongs to
            @showprogress for (i, element) in enumerate(lattice)
                if element.len > 0
                    # Calculate segments needed
                    segments_by_count = 10
                    segments_by_length = ceil(Int, element.len / 0.1)
                    num_segments = max(segments_by_count, segments_by_length)

                    # num_segments = 1

                    segment_length = element.len / num_segments

                    if isa(element, SBEND)
                        # Special handling for bends
                        angle_per_segment = element.angle / num_segments
                        
                        # First segment - keep entrance edge effect (e1)
                        first_segment = deepcopy(element)
                        first_segment.len = segment_length
                        first_segment.angle = angle_per_segment
                        if num_segments > 1
                            first_segment.e2 = 0.0  # No exit edge for first segment
                        end
                        push!(expanded_lattice, first_segment)
                        push!(expanded_indices, i)
                        
                        # Middle segments - no edge effects
                        for j in 2:num_segments-1
                            mid_segment = deepcopy(element)
                            mid_segment.len = segment_length
                            mid_segment.angle = angle_per_segment
                            mid_segment.e1 = 0.0
                            mid_segment.e2 = 0.0
                            push!(expanded_lattice, mid_segment)
                            push!(expanded_indices, i)
                        end
                        
                        # Last segment - keep exit edge effect (e2)
                        if num_segments > 1
                            last_segment = deepcopy(element)
                            last_segment.len = segment_length
                            last_segment.angle = angle_per_segment
                            last_segment.e1 = 0.0  # No entrance edge for last segment
                            push!(expanded_lattice, last_segment)
                            push!(expanded_indices, i)
                        end
                    else
                        # For other elements, proceed as before
                        for j in 1:num_segments
                            small_element = deepcopy(element)
                            small_element.len = segment_length
                            push!(expanded_lattice, small_element)
                            push!(expanded_indices, i)
                        end
                    end
                    
                    # # Create smaller elements
                    # for j in 1:num_segments
                    #     small_element = deepcopy(element)
                    #     small_element.len = segment_length
                    #     push!(expanded_lattice, small_element)
                    #     push!(expanded_indices, i)  # Record original index
                    # end

                    # push!(expanded_lattice, element)
                    # push!(expanded_indices, i)  # Record original index
                else
                    # Zero-length elements
                    push!(expanded_lattice, element)
                    push!(expanded_indices, i)
                end
            end
            
            # Propagate through expanded lattice
            n_expanded = length(expanded_lattice)
            expanded_beam_rms = zeros(n_expanded, 6)
            expanded_floor_distance = zeros(n_expanded)
            expanded_twi = zeros(n_expanded, 9)
            flat_particles = zeros(6 * np)
            mean_arr = zeros(n_expanded, 6)
            
            # Propagate through expanded lattice
            for i in eachindex(expanded_lattice)
                expanded_twi[i, :] .= twiss_beam(beam)
                mean_arr[i, :] .= vec(mean(beam.r, dims=1))
                flat_particles .= collect(Iterators.flatten(eachrow(beam.r)))
                
                # Pass particles through element
                pass!(expanded_lattice[i], flat_particles, np, beam)
                beam.r = reshape(flat_particles, 6, np)'
                
                # Store RMS values
                for dim in 1:6
                    expanded_beam_rms[i, dim] = sqrt(mean(beam.r[:,dim].^2))
                end
                
                # Calculate floor distance
                if i == 1
                    expanded_floor_distance[i] = expanded_lattice[i].len
                else
                    expanded_floor_distance[i] = expanded_floor_distance[i-1] + expanded_lattice[i].len
                end
            end
            
            # Return everything needed
            return expanded_beam_rms, expanded_floor_distance, beam, expanded_twi, 
                lattice, expanded_lattice, expanded_indices, mean_arr
        end

        function run_simulation()
            lattice, beam1, beam2, beam3 = create_lattice()
            
            beam1_rms, floor_distance, beam1, twi1, original_lattice, _, expanded_indices1 = 
                propagate_beam(lattice, beam1, beam1.np)
            
            beam2_rms, _, beam2, twi2, _, _, expanded_indices2 = 
                propagate_beam(lattice, beam2, beam2.np)
            
            beam3_rms, _, beam3, twi3, _, _, expanded_indices3 = 
                propagate_beam(lattice, beam3, beam3.np)

            return original_lattice, beam1, beam2, beam3, 
                beam1_rms, beam2_rms, beam3_rms, 
                twi1, twi2, twi3, floor_distance, expanded_indices1
        end

        # plotting wrapper function
        function plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, 
                                            twi1, twi2, twi3, floor_distance,
                                            original_lattice, expanded_indices, 
                                            combined=false)
            # Filter expanded data up to the specified original element
            mask = expanded_indices .<= end_ele
            
            # Use the filtered data for plots but original lattice for floor layout
            return plot_multibeam_data(
                floor_distance[mask], 
                beam1_rms[mask,:], 
                beam2_rms[mask,:], 
                beam3_rms[mask,:], 
                twi1[mask,:], 
                twi2[mask,:], 
                twi3[mask,:], 
                original_lattice[1:end_ele],  # Use original lattice for floor plot
                combined=combined
            )
        end

        function visualize_beam_properties(beam::Beam)
            # Create a multi-panel figure
            fig = Figure(size=(900, 600))
            
            # Phase space plots (x-px, y-py, z-dp)
            ax1 = Axis(fig[1, 1], title="x-px phase space", xlabel="x [mm]", ylabel="px")
            # Convert m to mm for x-axis
            scatter!(ax1, beam.r[:,1] .* 1000, beam.r[:,2], markersize=1, strokewidth=0, alpha=0.5)
            
            ax2 = Axis(fig[1, 2], title="y-py phase space", xlabel="y [mm]", ylabel="py")
            # Convert m to mm for y-axis
            scatter!(ax2, beam.r[:,3] .* 1000, beam.r[:,4], markersize=1, strokewidth=0, alpha=0.5)
            
            ax3 = Axis(fig[1, 3], title="z-dp phase space", xlabel="z [mm]", ylabel="dp/p")
            # Convert m to mm for z-axis
            scatter!(ax3, beam.r[:,5] .* 1000, beam.r[:,6], markersize=1, strokewidth=0, alpha=0.5)
            
            # Projections (histograms)
            ax4 = Axis(fig[2, 1], title="x distribution", xlabel="x [mm]")
            # Convert m to mm for x-axis
            hist!(ax4, beam.r[:,1] .* 1000, bins=50, color=(:blue, 0.7))
            
            ax5 = Axis(fig[2, 2], title="y distribution", xlabel="y [mm]")
            # Convert m to mm for y-axis
            hist!(ax5, beam.r[:,3] .* 1000, bins=50, color=(:blue, 0.7))
            
            ax6 = Axis(fig[2, 3], title="z distribution", xlabel="z [mm]")
            # Convert m to mm for z-axis
            hist!(ax6, beam.r[:,5] .* 1000, bins=50, color=(:blue, 0.7))
            
            # Return the figure
            return fig
        end

        # function plot_multibeam_data(floor_length, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, lattice; combined = true)
        #     if combined != true    
        #         # Define beam colors and labels (shared across all plots)
        #         beam_colors = [:blue, :orange, :green]
        #         beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]
                
        #         # Define element colors and heights (shared across all plots)
        #         element_colors = Dict(
        #             "SBEND" => :purple,
        #             "KQUAD" => :red,
        #             "ORBTRIM" => :green,
        #             "MARKER" => :gray,
        #             "DRIFT" => :white,
        #             "default" => :gray
        #         )
                
        #         element_heights = Dict(
        #             "SBEND" => 0.8,
        #             "KQUAD" => 0.7,
        #             "ORBTRIM" => 1.0,
        #             "MARKER" => 0.4,
        #             "DRIFT" => 0.3,
        #             "default" => 0.5
        #         )
                
        #         # Helper function to create a standard plot structure
        #         function create_plot(top_data1, top_data2, top_data3, bottom_data1, bottom_data2, bottom_data3, 
        #                             top_label, bottom_label, title)
                    
        #             f = Figure(size = (1200, 400))
        #             gl = f[1, 1] = GridLayout()
        #             legend_layout = f[1, 2] = GridLayout()
                    
        #             # Top plot
        #             ax_top = Axis(gl[1, 1], xlabel = "", ylabel = top_label)
        #             lines!(ax_top, floor_length, top_data1, color = beam_colors[1], linewidth = 2)
        #             lines!(ax_top, floor_length, top_data2, color = beam_colors[2], linewidth = 2)
        #             lines!(ax_top, floor_length, top_data3, color = beam_colors[3], linewidth = 2)
        #             hidexdecorations!(ax_top)
                    
        #             # Create floorline plot in the middle
        #             floor_ax = Axis(gl[2, 1], 
        #                         xlabel = "",
        #                         ylabel = "",
        #                         yticklabelsvisible = false,
        #                         yticksvisible = false)
                    
        #             # Track element types for legend
        #             element_types = Dict{String, Any}()
                    
        #             # Draw the floorline
        #             curr_pos = 0.0
        #             for ele in lattice
        #                 # Get element type
        #                 ele_type = string(typeof(ele).name.name)
                        
        #                 # Determine color and height
        #                 color = get(element_colors, ele_type, element_colors["default"])
        #                 height = get(element_heights, ele_type, element_heights["default"])
                        
        #                 # Make sure elements don't overlap by using a minimum width
        #                 ele_width = ele.len #max(ele.len, 0.05)
                        
        #                 # Draw rectangle for element based on element type and properties
        #                 if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
        #                     # Use the original quad height
        #                     quad_height = element_heights["KQUAD"]
        #                     drift_height = element_heights["DRIFT"]
                            
        #                     # Calculate offset to align with other elements
        #                     # This will make it overlap the x-axis slightly to maintain alignment
        #                     if ele.k1 > 0
        #                         # For focusing quads: align bottom with the bottom of drift elements
        #                         y_pos = -drift_height/2
        #                     else
        #                         # For defocusing quads: align top with the top of drift elements
        #                         y_pos = drift_height/2 - quad_height
        #                     end
                            
        #                     rect = Rect(curr_pos, y_pos, ele_width, quad_height)
        #                 elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
        #                     # Use the original sextupole height
        #                     sext_height = get(element_heights, "KSEXT", element_heights["default"])
        #                     drift_height = element_heights["DRIFT"]
                            
        #                     # Calculate offset to align with other elements
        #                     if ele.k2 > 0
        #                         # For focusing sextupoles: align bottom with the bottom of drift elements
        #                         y_pos = -drift_height/2
        #                     else
        #                         # For defocusing sextupoles: align top with the top of drift elements
        #                         y_pos = drift_height/2 - sext_height
        #                     end
                            
        #                     rect = Rect(curr_pos, y_pos, ele_width, sext_height)
        #                 else
        #                     # All other elements centered on x-axis as before
        #                     rect = Rect(curr_pos, -height/2, ele_width, height)
        #                 end
                        
        #                 element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                        
        #                 # Track for legend (only if not already tracked)
        #                 if !haskey(element_types, ele_type)
        #                     element_types[ele_type] = element
        #                 end
                        
        #                 # Update position - ensure we advance by at least the element width
        #                 # to avoid overlapping due to position calculations
        #                 curr_pos += ele_width
        #             end
                    
        #             hidexdecorations!(floor_ax)
                    
        #             # Bottom plot
        #             ax_bottom = Axis(gl[3, 1], xlabel = "z [m]", ylabel = bottom_label, yreversed = true)
        #             lines!(ax_bottom, floor_length, bottom_data1, color = beam_colors[1], linewidth = 2)
        #             lines!(ax_bottom, floor_length, bottom_data2, color = beam_colors[2], linewidth = 2)
        #             lines!(ax_bottom, floor_length, bottom_data3, color = beam_colors[3], linewidth = 2)
                    
        #             # Link x axes
        #             linkxaxes!(ax_top, floor_ax, ax_bottom)
                    
        #             # Add beam legend
        #             legend_entries = [
        #                 LineElement(color = c, linewidth = 2) for c in beam_colors
        #             ]
                    
        #             leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                    
        #             # Add element legend (filtering out MARKER and DRIFT)
        #             element_entries = []
        #             element_names = []
        #             sorted_names = sort(collect(keys(element_types)))
        #             for name in sorted_names
        #                 if !(name in ["MARKER", "DRIFT"])
        #                     element = element_types[name]
        #                     push!(element_entries, element)
        #                     push!(element_names, name)
        #                 end
        #             end
                    
        #             if !isempty(element_entries)
        #                 element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
        #             end
                    
        #             # Adjust spacing between legends
        #             rowgap!(legend_layout, 5)
                    
        #             # Set row sizes
        #             rowsize!(gl, 1, 100)  # Top plot
        #             rowsize!(gl, 2, 30)   # Elements plot (smaller)
        #             rowsize!(gl, 3, 100)  # Bottom plot
                    
        #             # Set title
        #             # f.title = title
                    
        #             return f
        #         end
                
        #         # Create the four plots
        #         rms_plot = create_plot(
        #             beam1_rms[:,1] .* 1e3,
        #             beam2_rms[:,1] .* 1e3,
        #             beam3_rms[:,1] .* 1e3,
        #             beam1_rms[:,3] .* 1e3,
        #             beam2_rms[:,3] .* 1e3,
        #             beam3_rms[:,3] .* 1e3,
        #             "RMS X [mm]",
        #             "RMS Y [mm]",
        #             "RMS Plot"
        #         )
                
        #         beta_plot = create_plot(
        #             twi1[:,1],
        #             twi2[:,1],
        #             twi3[:,1],
        #             twi1[:,4],
        #             twi2[:,4],
        #             twi3[:,4],
        #             L"\beta_x",
        #             L"\beta_y",
        #             "Beta Functions"
        #         )
                
        #         alpha_plot = create_plot(
        #             twi1[:,2],
        #             twi2[:,2],
        #             twi3[:,2],
        #             twi1[:,5],
        #             twi2[:,5],
        #             twi3[:,5],
        #             L"\alpha_x",
        #             L"\alpha_y",
        #             "Alpha Functions"
        #         )
                
        #         emittance_plot = create_plot(
        #             twi1[:,3],
        #             twi2[:,3],
        #             twi3[:,3],
        #             twi1[:,6],
        #             twi2[:,6],
        #             twi3[:,6],
        #             L"\epsilon_x",
        #             L"\epsilon_y",
        #             "Emittance"
        #         )
                
        #         # Return all plots
        #         return Dict(
        #             "rms" => rms_plot,
        #             "beta" => beta_plot,
        #             "alpha" => alpha_plot,
        #             "emittance" => emittance_plot
        #         )
        #     else
        #         # Create one large figure
        #         f = Figure(size = (1200, 1200))
                
        #         # Create the main layout
        #         gl = f[1, 1] = GridLayout()
        #         legend_layout = f[1, 2] = GridLayout()
                
        #         # Define beam colors and labels
        #         beam_colors = [:blue, :orange, :green]
        #         beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]
                
        #         # Define element colors and heights
        #         element_colors = Dict(
        #             "SBEND" => :purple,
        #             "KQUAD" => :red,
        #             "ORBTRIM" => :green,
        #             "MARKER" => :gray,
        #             "DRIFT" => :white,
        #             "default" => :gray
        #         )
                
        #         element_heights = Dict(
        #             "SBEND" => 0.8,
        #             "KQUAD" => 0.7,
        #             "ORBTRIM" => 0.6,
        #             "MARKER" => 0.4,
        #             "DRIFT" => 0.3,
        #             "default" => 0.5
        #         )
                
        #         # Create plots in sequence (8 total + 1 lattice in the middle)
        #         # 1. RMS X
        #         ax_rms_x = Axis(gl[1, 1], xlabel = "", ylabel = "RMS X [mm]", title = "RMS Values")
        #         lines!(ax_rms_x, floor_length, beam1_rms[:,1] .* 1e3, color = beam_colors[1], linewidth = 2)
        #         lines!(ax_rms_x, floor_length, beam2_rms[:,1] .* 1e3, color = beam_colors[2], linewidth = 2)
        #         lines!(ax_rms_x, floor_length, beam3_rms[:,1] .* 1e3, color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_rms_x)
                
        #         # 2. RMS Y
        #         ax_rms_y = Axis(gl[2, 1], xlabel = "", ylabel = "RMS Y [mm]")
        #         lines!(ax_rms_y, floor_length, beam1_rms[:,3] .* 1e3, color = beam_colors[1], linewidth = 2)
        #         lines!(ax_rms_y, floor_length, beam2_rms[:,3] .* 1e3, color = beam_colors[2], linewidth = 2)
        #         lines!(ax_rms_y, floor_length, beam3_rms[:,3] .* 1e3, color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_rms_y)
                
        #         # 3. Beta X
        #         ax_beta_x = Axis(gl[3, 1], xlabel = "", ylabel = L"\beta_x", title = "Beta Functions")
        #         lines!(ax_beta_x, floor_length, twi1[:,1], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_beta_x, floor_length, twi2[:,1], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_beta_x, floor_length, twi3[:,1], color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_beta_x)
                
        #         # 4. Beta Y
        #         ax_beta_y = Axis(gl[4, 1], xlabel = "", ylabel = L"\beta_y")
        #         lines!(ax_beta_y, floor_length, twi1[:,4], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_beta_y, floor_length, twi2[:,4], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_beta_y, floor_length, twi3[:,4], color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_beta_y)
                
        #         # 5. Lattice plot in the middle
        #         floor_ax = Axis(gl[5, 1], 
        #                     xlabel = "",
        #                     ylabel = "",
        #                     yticklabelsvisible = false,
        #                     yticksvisible = false)
                
        #         # Track element types for legend
        #         element_types = Dict{String, Any}()
                
        #         # Draw the floorline
        #         curr_pos = 0.0
        #         for ele in lattice
        #             # Get element type
        #             ele_type = string(typeof(ele).name.name)
                    
        #             # Determine color and height
        #             color = get(element_colors, ele_type, element_colors["default"])
        #             height = get(element_heights, ele_type, element_heights["default"])
                    
        #             # Use a minimum width to prevent tiny elements and avoid overlap
        #             ele_width = ele.len #max(ele.len, 0.05)
                    
        #             # Draw rectangle for element based on element type and properties
        #             if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
        #                 # Use the original quad height, not drift height
        #                 quad_height = element_heights["KQUAD"]
        #                 if ele.k1 > 0
        #                     # Position touching x-axis and extending upward with full quad height
        #                     rect = Rect(curr_pos, 0, ele_width, quad_height)
        #                 else
        #                     # Position touching x-axis and extending downward with full quad height
        #                     rect = Rect(curr_pos, -quad_height, ele_width, quad_height)
        #                 end
        #             elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
        #                 # Use the original sextupole height
        #                 sext_height = get(element_heights, "KSEXT", element_heights["default"])
        #                 if ele.k2 > 0
        #                     # Position touching x-axis and extending upward with full sextupole height
        #                     rect = Rect(curr_pos, 0, ele_width, sext_height)
        #                 else
        #                     # Position touching x-axis and extending downward with full sextupole height
        #                     rect = Rect(curr_pos, -sext_height, ele_width, sext_height)
        #                 end
        #             else
        #                 # All other elements centered on x-axis as before
        #                 rect = Rect(curr_pos, -height/2, ele_width, height)
        #             end
                    
        #             element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                    
        #             # Track for legend (only if not already tracked)
        #             if !haskey(element_types, ele_type)
        #                 element_types[ele_type] = element
        #             end
                    
        #             # Update position - use the element width to ensure no overlap
        #             curr_pos += ele_width
        #         end
        #         hidexdecorations!(floor_ax)
                
        #         # 6. Alpha X
        #         ax_alpha_x = Axis(gl[6, 1], xlabel = "", ylabel = L"\alpha_x", title = "Alpha Functions")
        #         lines!(ax_alpha_x, floor_length, twi1[:,2], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_alpha_x, floor_length, twi2[:,2], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_alpha_x, floor_length, twi3[:,2], color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_alpha_x)
                
        #         # 7. Alpha Y
        #         ax_alpha_y = Axis(gl[7, 1], xlabel = "", ylabel = L"\alpha_y")
        #         lines!(ax_alpha_y, floor_length, twi1[:,5], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_alpha_y, floor_length, twi2[:,5], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_alpha_y, floor_length, twi3[:,5], color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_alpha_y)
                
        #         # 8. Emittance X
        #         ax_emit_x = Axis(gl[8, 1], xlabel = "", ylabel = L"\epsilon_x", title = "Emittance")
        #         lines!(ax_emit_x, floor_length, twi1[:,3], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_emit_x, floor_length, twi2[:,3], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_emit_x, floor_length, twi3[:,3], color = beam_colors[3], linewidth = 2)
        #         hidexdecorations!(ax_emit_x)
                
        #         # 9. Emittance Y (bottom plot)
        #         ax_emit_y = Axis(gl[9, 1], xlabel = "z [m]", ylabel = L"\epsilon_y")
        #         lines!(ax_emit_y, floor_length, twi1[:,6], color = beam_colors[1], linewidth = 2)
        #         lines!(ax_emit_y, floor_length, twi2[:,6], color = beam_colors[2], linewidth = 2)
        #         lines!(ax_emit_y, floor_length, twi3[:,6], color = beam_colors[3], linewidth = 2)
                
        #         # Link all x axes
        #         for ax in [ax_rms_x, ax_rms_y, ax_beta_x, ax_beta_y, floor_ax, 
        #                 ax_alpha_x, ax_alpha_y, ax_emit_x, ax_emit_y]
        #             linkxaxes!(ax, ax_rms_x)
        #         end
                
        #         # Add beam legend
        #         legend_entries = [
        #             LineElement(color = c, linewidth = 2) for c in beam_colors
        #         ]
                
        #         leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                
        #         # Add element legend (excluding MARKER and DRIFT)
        #         element_entries = []
        #         element_names = []
        #         sorted_names = sort(collect(keys(element_types)))
        #         for name in sorted_names
        #             if !(name in ["MARKER", "DRIFT"])
        #                 element = element_types[name]
        #                 push!(element_entries, element)
        #                 push!(element_names, name)
        #             end
        #         end
                
        #         if !isempty(element_entries)
        #             element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
        #         end
                
        #         # Adjust spacing between legends
        #         rowgap!(legend_layout, 5)
                
        #         # Set row sizes (make lattice row smaller)
        #         rowsize!(gl, 5, 30)  # Lattice plot is smaller
                
                
        #         return f
        #     end
        # end

        function plot_multibeam_data(floor_length, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, lattice; combined = true)
            if combined != true    
                # Define beam colors and labels (shared across all plots)
                # beam_colors = [:blue , :orange, :green]
                # beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]

                beam_colors = [:blue]
                beam_labels = [L"^{124}Xe^{50+}"]
                
                # Define element colors and heights (shared across all plots)
                element_colors = Dict(
                    "SBEND" => :purple,
                    "KQUAD" => :red,
                    "ORBTRIM" => :green,
                    "MARKER" => :gray,
                    "DRIFT" => :white,
                    "default" => :gray
                )
                
                element_heights = Dict(
                    "SBEND" => 0.8,
                    "KQUAD" => 0.7,
                    "ORBTRIM" => 1.0,
                    "MARKER" => 0.4,
                    "DRIFT" => 0.3,
                    "default" => 0.5
                )
                
                # Helper function to create a standard plot structure
                function create_plot(top_data1, top_data2, top_data3, bottom_data1, bottom_data2, bottom_data3, 
                                    top_label, bottom_label, title)
                    
                    f = Figure(size = (1200, 400))
                    gl = f[1, 1] = GridLayout()
                    legend_layout = f[1, 2] = GridLayout()
                    
                    # Top plot
                    ax_top = Axis(gl[1, 1], xlabel = "", ylabel = top_label)
                    lines!(ax_top, floor_length, top_data1, color = beam_colors[1], linewidth = 2)
                    # lines!(ax_top, floor_length, top_data2, color = beam_colors[2], linewidth = 2)
                    # lines!(ax_top, floor_length, top_data3, color = beam_colors[3], linewidth = 2)
                    hidexdecorations!(ax_top)
                    
                    # Create floorline plot in the middle
                    floor_ax = Axis(gl[2, 1], 
                                xlabel = "",
                                ylabel = "",
                                yticklabelsvisible = false,
                                yticksvisible = false)
                    
                    # Track element types for legend
                    element_types = Dict{String, Any}()
                    
                    # Draw the floorline
                    curr_pos = 0.0
                    for ele in lattice
                        # Get element type
                        ele_type = string(typeof(ele).name.name)
                        
                        # Determine color and height
                        color = get(element_colors, ele_type, element_colors["default"])
                        height = get(element_heights, ele_type, element_heights["default"])
                        
                        # Make sure elements don't overlap by using a minimum width
                        ele_width = ele.len #max(ele.len, 0.05)
                        
                        # Draw rectangle for element based on element type and properties
                        if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
                            # Use the original quad height
                            quad_height = element_heights["KQUAD"]
                            drift_height = element_heights["DRIFT"]
                            
                            # Calculate offset to align with other elements
                            # This will make it overlap the x-axis slightly to maintain alignment
                            if ele.k1 > 0
                                # For focusing quads: align bottom with the bottom of drift elements
                                y_pos = -drift_height/2
                            else
                                # For defocusing quads: align top with the top of drift elements
                                y_pos = drift_height/2 - quad_height
                            end
                            
                            rect = Rect(curr_pos, y_pos, ele_width, quad_height)
                        elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
                            # Use the original sextupole height
                            sext_height = get(element_heights, "KSEXT", element_heights["default"])
                            drift_height = element_heights["DRIFT"]
                            
                            # Calculate offset to align with other elements
                            if ele.k2 > 0
                                # For focusing sextupoles: align bottom with the bottom of drift elements
                                y_pos = -drift_height/2
                            else
                                # For defocusing sextupoles: align top with the top of drift elements
                                y_pos = drift_height/2 - sext_height
                            end
                            
                            rect = Rect(curr_pos, y_pos, ele_width, sext_height)
                        else
                            # All other elements centered on x-axis as before
                            rect = Rect(curr_pos, -height/2, ele_width, height)
                        end
                        
                        element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                        
                        # Track for legend (only if not already tracked)
                        if !haskey(element_types, ele_type)
                            element_types[ele_type] = element
                        end
                        
                        # Update position - ensure we advance by at least the element width
                        # to avoid overlapping due to position calculations
                        curr_pos += ele_width
                    end
                    
                    hidexdecorations!(floor_ax)
                    
                    # Bottom plot
                    ax_bottom = Axis(gl[3, 1], xlabel = "z [m]", ylabel = bottom_label, yreversed = true)
                    lines!(ax_bottom, floor_length, bottom_data1, color = beam_colors[1], linewidth = 2)
                    # lines!(ax_bottom, floor_length, bottom_data2, color = beam_colors[2], linewidth = 2)
                    # lines!(ax_bottom, floor_length, bottom_data3, color = beam_colors[3], linewidth = 2)
                    
                    # Link x axes
                    linkxaxes!(ax_top, floor_ax, ax_bottom)
                    
                    # Add beam legend
                    legend_entries = [
                        LineElement(color = c, linewidth = 2) for c in beam_colors
                    ]
                    
                    leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                    
                    # Add element legend (filtering out MARKER and DRIFT)
                    element_entries = []
                    element_names = []
                    sorted_names = sort(collect(keys(element_types)))
                    for name in sorted_names
                        if !(name in ["MARKER", "DRIFT"])
                            element = element_types[name]
                            push!(element_entries, element)
                            push!(element_names, name)
                        end
                    end
                    
                    if !isempty(element_entries)
                        element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
                    end
                    
                    # Adjust spacing between legends
                    rowgap!(legend_layout, 5)
                    
                    # Set row sizes
                    rowsize!(gl, 1, 100)  # Top plot
                    rowsize!(gl, 2, 30)   # Elements plot (smaller)
                    rowsize!(gl, 3, 100)  # Bottom plot
                    
                    # Set title
                    # f.title = title
                    
                    return f
                end
                
                # Create the four plots
                rms_plot = create_plot(
                    beam1_rms[:,1] .* 1e3,
                    beam2_rms[:,1] .* 1e3,
                    beam3_rms[:,1] .* 1e3,
                    beam1_rms[:,3] .* 1e3,
                    beam2_rms[:,3] .* 1e3,
                    beam3_rms[:,3] .* 1e3,
                    "RMS X [mm]",
                    "RMS Y [mm]",
                    "RMS Plot"
                )
                
                beta_plot = create_plot(
                    twi1[:,1],
                    twi2[:,1],
                    twi3[:,1],
                    twi1[:,4],
                    twi2[:,4],
                    twi3[:,4],
                    L"\beta_x",
                    L"\beta_y",
                    "Beta Functions"
                )
                
                alpha_plot = create_plot(
                    twi1[:,2],
                    twi2[:,2],
                    twi3[:,2],
                    twi1[:,5],
                    twi2[:,5],
                    twi3[:,5],
                    L"\alpha_x",
                    L"\alpha_y",
                    "Alpha Functions"
                )
                
                emittance_plot = create_plot(
                    twi1[:,3],
                    twi2[:,3],
                    twi3[:,3],
                    twi1[:,6],
                    twi2[:,6],
                    twi3[:,6],
                    L"\epsilon_x",
                    L"\epsilon_y",
                    "Emittance"
                )
                
                # Return all plots
                return Dict(
                    "rms" => rms_plot,
                    "beta" => beta_plot,
                    "alpha" => alpha_plot,
                    "emittance" => emittance_plot
                )
            else
                # Create one large figure
                f = Figure(size = (1200, 1200))
                
                # Create the main layout
                gl = f[1, 1] = GridLayout()
                legend_layout = f[1, 2] = GridLayout()
                
                # Define beam colors and labels
                # beam_colors = [:blue , :orange, :green]
                # beam_labels = [L"^{124}Xe^{50+}", L"^{124}Xe^{49+}", L"^{124}Xe^{51+}"]

                beam_colors = [:blue]
                beam_labels = [L"^{124}Xe^{50+}"]
                
                # Define element colors and heights
                element_colors = Dict(
                    "SBEND" => :purple,
                    "KQUAD" => :red,
                    "ORBTRIM" => :green,
                    "MARKER" => :gray,
                    "DRIFT" => :white,
                    "default" => :gray
                )
                
                element_heights = Dict(
                    "SBEND" => 0.8,
                    "KQUAD" => 0.7,
                    "ORBTRIM" => 0.6,
                    "MARKER" => 0.4,
                    "DRIFT" => 0.3,
                    "default" => 0.5
                )
                
                # Create plots in sequence (8 total + 1 lattice in the middle)
                # 1. RMS X
                ax_rms_x = Axis(gl[1, 1], xlabel = "", ylabel = "RMS X [mm]", title = "RMS Values")
                lines!(ax_rms_x, floor_length, beam1_rms[:,1] .* 1e3, color = beam_colors[1], linewidth = 2)
                # lines!(ax_rms_x, floor_length, beam2_rms[:,1] .* 1e3, color = beam_colors[2], linewidth = 2)
                # lines!(ax_rms_x, floor_length, beam3_rms[:,1] .* 1e3, color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_rms_x)
                
                # 2. RMS Y
                ax_rms_y = Axis(gl[2, 1], xlabel = "", ylabel = "RMS Y [mm]")
                lines!(ax_rms_y, floor_length, beam1_rms[:,3] .* 1e3, color = beam_colors[1], linewidth = 2)
                # lines!(ax_rms_y, floor_length, beam2_rms[:,3] .* 1e3, color = beam_colors[2], linewidth = 2)
                # lines!(ax_rms_y, floor_length, beam3_rms[:,3] .* 1e3, color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_rms_y)
                
                # 3. Beta X
                ax_beta_x = Axis(gl[3, 1], xlabel = "", ylabel = L"\beta_x", title = "Beta Functions")
                lines!(ax_beta_x, floor_length, twi1[:,1], color = beam_colors[1], linewidth = 2)
                # lines!(ax_beta_x, floor_length, twi2[:,1], color = beam_colors[2], linewidth = 2)
                # lines!(ax_beta_x, floor_length, twi3[:,1], color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_beta_x)
                
                # 4. Beta Y
                ax_beta_y = Axis(gl[4, 1], xlabel = "", ylabel = L"\beta_y")
                lines!(ax_beta_y, floor_length, twi1[:,4], color = beam_colors[1], linewidth = 2)
                # lines!(ax_beta_y, floor_length, twi2[:,4], color = beam_colors[2], linewidth = 2)
                # lines!(ax_beta_y, floor_length, twi3[:,4], color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_beta_y)
                
                # 5. Lattice plot in the middle
                floor_ax = Axis(gl[5, 1], 
                            xlabel = "",
                            ylabel = "",
                            yticklabelsvisible = false,
                            yticksvisible = false)
                
                # Track element types for legend
                element_types = Dict{String, Any}()
                
                # Draw the floorline
                curr_pos = 0.0
                for ele in lattice
                    # Get element type
                    ele_type = string(typeof(ele).name.name)
                    
                    # Determine color and height
                    color = get(element_colors, ele_type, element_colors["default"])
                    height = get(element_heights, ele_type, element_heights["default"])
                    
                    # Use a minimum width to prevent tiny elements and avoid overlap
                    ele_width = ele.len #max(ele.len, 0.05)
                    
                    # Draw rectangle for element based on element type and properties
                    if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
                        # Use the original quad height, not drift height
                        quad_height = element_heights["KQUAD"]
                        if ele.k1 > 0
                            # Position touching x-axis and extending upward with full quad height
                            rect = Rect(curr_pos, 0, ele_width, quad_height)
                        else
                            # Position touching x-axis and extending downward with full quad height
                            rect = Rect(curr_pos, -quad_height, ele_width, quad_height)
                        end
                    elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
                        # Use the original sextupole height
                        sext_height = get(element_heights, "KSEXT", element_heights["default"])
                        if ele.k2 > 0
                            # Position touching x-axis and extending upward with full sextupole height
                            rect = Rect(curr_pos, 0, ele_width, sext_height)
                        else
                            # Position touching x-axis and extending downward with full sextupole height
                            rect = Rect(curr_pos, -sext_height, ele_width, sext_height)
                        end
                    else
                        # All other elements centered on x-axis as before
                        rect = Rect(curr_pos, -height/2, ele_width, height)
                    end
                    
                    element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                    
                    # Track for legend (only if not already tracked)
                    if !haskey(element_types, ele_type)
                        element_types[ele_type] = element
                    end
                    
                    # Update position - use the element width to ensure no overlap
                    curr_pos += ele_width
                end
                hidexdecorations!(floor_ax)
                
                # 6. Alpha X
                ax_alpha_x = Axis(gl[6, 1], xlabel = "", ylabel = L"\alpha_x", title = "Alpha Functions")
                lines!(ax_alpha_x, floor_length, twi1[:,2], color = beam_colors[1], linewidth = 2)
                # lines!(ax_alpha_x, floor_length, twi2[:,2], color = beam_colors[2], linewidth = 2)
                # lines!(ax_alpha_x, floor_length, twi3[:,2], color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_alpha_x)
                
                # 7. Alpha Y
                ax_alpha_y = Axis(gl[7, 1], xlabel = "", ylabel = L"\alpha_y")
                lines!(ax_alpha_y, floor_length, twi1[:,5], color = beam_colors[1], linewidth = 2)
                # lines!(ax_alpha_y, floor_length, twi2[:,5], color = beam_colors[2], linewidth = 2)
                # lines!(ax_alpha_y, floor_length, twi3[:,5], color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_alpha_y)
                
                # 8. Emittance X
                ax_emit_x = Axis(gl[8, 1], xlabel = "", ylabel = L"\epsilon_x", title = "Emittance")
                lines!(ax_emit_x, floor_length, twi1[:,3], color = beam_colors[1], linewidth = 2)
                # lines!(ax_emit_x, floor_length, twi2[:,3], color = beam_colors[2], linewidth = 2)
                # lines!(ax_emit_x, floor_length, twi3[:,3], color = beam_colors[3], linewidth = 2)
                hidexdecorations!(ax_emit_x)
                
                # 9. Emittance Y (bottom plot)
                ax_emit_y = Axis(gl[9, 1], xlabel = "z [m]", ylabel = L"\epsilon_y")
                lines!(ax_emit_y, floor_length, twi1[:,6], color = beam_colors[1], linewidth = 2)
                # lines!(ax_emit_y, floor_length, twi2[:,6], color = beam_colors[2], linewidth = 2)
                # lines!(ax_emit_y, floor_length, twi3[:,6], color = beam_colors[3], linewidth = 2)
                
                # Link all x axes
                for ax in [ax_rms_x, ax_rms_y, ax_beta_x, ax_beta_y, floor_ax, 
                        ax_alpha_x, ax_alpha_y, ax_emit_x, ax_emit_y]
                    linkxaxes!(ax, ax_rms_x)
                end
                
                # Add beam legend
                legend_entries = [
                    LineElement(color = c, linewidth = 2) for c in beam_colors
                ]
                
                leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
                
                # Add element legend (excluding MARKER and DRIFT)
                element_entries = []
                element_names = []
                sorted_names = sort(collect(keys(element_types)))
                for name in sorted_names
                    if !(name in ["MARKER", "DRIFT"])
                        element = element_types[name]
                        push!(element_entries, element)
                        push!(element_names, name)
                    end
                end
                
                if !isempty(element_entries)
                    element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
                end
                
                # Adjust spacing between legends
                rowgap!(legend_layout, 5)
                
                # Set row sizes (make lattice row smaller)
                rowsize!(gl, 5, 30)  # Lattice plot is smaller
                
                
                return f
            end
        end

        function reconstruct_original_floor_distance(expanded, distance)
            # Find the maximum original index
            max_index = maximum(expanded)
            
            # Initialize the original vector with zeros
            original = zeros(eltype(distance), max_index)
            

            for i in 1:length(expanded)
                original[expanded[i]] = distance[i]
            end
            
            return original
        end

        function FLAME_to_Julia_matrix(matrix)
            new_mat = reshape(matrix, 7, 7)
            return new_mat[1:end-1, 1:end-1]
        end

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
        
        # Create the lattice and beams
        lattice, beam1, beam2, beam3 = create_lattice();
        

        begin
            p1 = visualize_beam_properties(beam1);
            p2 = visualize_beam_properties(beam2);
            p3 = visualize_beam_properties(beam3);
            display(p1)
            display(p2)
            display(p3)
        end;


        # Run simulation
        lattice, beam1, beam2, beam3, beam1_rms, beam2_rms, beam3_rms, 
        twi1, twi2, twi3, floor_distance, expanded_indices = run_simulation();

        # Plot using original indexing
        end_ele = 200;
        plots_jutrack = plot_multibeam(end_ele, beam1_centered_rms, beam2_centered_rms, beam3_centered_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, false);
        plots_jutrack["rms"]
        save("src/demo/FRIB/plots/rms_first_50.png", plots_jutrack["rms"])
        p4 = plot_multibeam(end_ele, beam1_centered_rms, beam2_centered_rms, beam3_centered_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, true)
        save("src/demo/FRIB/plots/all_first_50.png", p4)

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

    # if mapped_name == "k1"
    #     param_value = param_value / 5.69860938
    # end
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

