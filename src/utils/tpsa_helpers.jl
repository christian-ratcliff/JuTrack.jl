using YAML

function load_flame_data(yaml_file::String; charges::Vector{Int}=[49, 50, 51])
    yaml_data = YAML.load_file(yaml_file)
    lattice_data = yaml_data["lattice"]
    settings = lattice_data["settings"]
    beams_data = lattice_data["beams"]
    
    # Extract settings for coordinate transformation
    rest_energy = settings["IonEs"]
    kinetic_energy = settings["IonEk"]
    mass_number = settings["mass_number"]
    rf_freq = settings["rf_freq"]
    brho_values = Dict{Int, Float64}()
    for (k, v) in settings["brho_values"]
        brho_values[parse(Int, string(k))] = Float64(v)
    end
    reference_charge = settings["reference_charge"]
    
    # Calculate transformation parameters
    gamma = 1.0 + kinetic_energy / rest_energy
    beta = sqrt(1.0 - 1.0 / gamma^2)
    rf_k = 2 * pi * rf_freq / 299792458.0
    scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(rest_energy+kinetic_energy)]
    
    # Extract and transform beam data for each charge
    flame_matrices = Dict{Int, Matrix{Float64}}()
    flame_centroids = Dict{Int, Vector{Float64}}()
    
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
            error("No beam data found for charge $charge")
        end
        
        # Get raw FLAME data
        matrix_flat = beam_data["matrix"]
        matrix = reshape(matrix_flat, 6, 6)
        centroid = beam_data["centroid"]
        
        # Transform FLAME matrix to JuTrack matrix
        matrix_transformed = diagm(scaling) * ((matrix + matrix') ./ 2) * diagm(scaling)'
        
        # Transform centroid with scaling
        centroid_transformed = centroid .* scaling
        
        # Calculate and add momentum deviation
        reference_brho = brho_values[reference_charge]
        beam_brho = brho_values[charge]
        momentum_deviation = (beam_brho - reference_brho) / reference_brho
        centroid_transformed[6] += momentum_deviation
        
        flame_matrices[charge] = matrix_transformed
        flame_centroids[charge] = centroid_transformed
    end
    
    return flame_matrices, flame_centroids, settings
end

function create_tpsa_beam(centroid)
    tpsa_coords = Vector{CTPS{Float64, 6, 1}}(undef, 6)
    for i in 1:6
        tpsa_coords[i] = CTPS(centroid[i], i, 6, 1)
    end
    return tpsa_coords
end

function extract_transfer_matrix(tpsa_coords)
    matrix = zeros(6, 6)
    for i in 1:6
        for j in 1:6
            matrix[i, j] = tpsa_coords[i].map[j+1]  
        end
    end
    return matrix
end

function extract_centroid(tpsa_coords)
    centroid = zeros(6)
    for i in 1:6
        centroid[i] = tpsa_coords[i].map[1]  
    end
    return centroid
end