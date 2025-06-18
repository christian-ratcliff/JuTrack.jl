begin
    include("../../JuTrack.jl")
    using .JuTrack
    using Statistics
    using LinearAlgebra
    using DelimitedFiles
    using YAML
    using BenchmarkTools
end

begin
    function benchmark_tpsa_method()
        # Load raw FLAME data from YAML
        flame_matrices, flame_centroids, settings = load_flame_data(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        # Load lattice
        lattice, _ = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        charges = sort(collect(keys(flame_matrices)))
        
        # Propagate FLAME data using TPSA
        transfer_matrices, beam_matrices, centroids_tpsa, floor_distance = 
            propagate_flame_data_TPSA(lattice, flame_matrices, flame_centroids, charges)
        
        # Extract beam properties from matrices
        beam_rms, beam_mean, beam_centered_rms, twi = 
            extract_beam_properties_from_matrices(beam_matrices, centroids_tpsa, floor_distance, charges)
        
        return beam_rms, beam_mean, beam_centered_rms, twi, floor_distance
    end

    function benchmark_particle_method()
        # Load lattice and create particle distributions from FLAME data
        lattice, multibeam = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        n_elements = length(lattice)
        charges = sort(collect(keys(multibeam.beams)))
        
        beam_rms = Dict{Int, Matrix{Float64}}()
        beam_mean = Dict{Int, Matrix{Float64}}()
        beam_centered_rms = Dict{Int, Matrix{Float64}}()
        twi = Dict{Int, Matrix{Float64}}()
        
        for charge in charges
            beam_rms[charge] = zeros(n_elements, 6)
            beam_mean[charge] = zeros(n_elements, 6)
            beam_centered_rms[charge] = zeros(n_elements, 6)
            twi[charge] = zeros(n_elements, 9)
        end
        
        floor_distance = zeros(n_elements)
        
        # Track particles through lattice
        for i in eachindex(lattice)
            pass!(lattice[i], multibeam)
            
            for charge in charges
                beam = multibeam.beams[charge]
                
                twi[charge][i, :] .= twiss_beam(beam)
                means = mean(beam.r, dims=1)
                beam_mean[charge][i, :] .= vec(means)
                
                for dim in 1:6
                    beam_rms[charge][i, dim] = sqrt(mean(beam.r[:,dim].^2))
                    beam_centered_rms[charge][i, dim] = sqrt(mean((beam.r[:,dim] .- means[dim]).^2))
                end
            end
            
            if i == 1
                floor_distance[i] = lattice[i].len
            else
                floor_distance[i] = floor_distance[i-1] + lattice[i].len
            end
        end
        
        return beam_rms, beam_mean, beam_centered_rms, twi, floor_distance
    end

    function run_performance_comparison()
        
        println("=== Performance Benchmark Comparison ===")
        println("Running benchmarks for full lattice simulation...")
        
        # Benchmark TPSA method
        println("Benchmarking TPSA method...")
        tpsa_benchmark = @benchmark benchmark_tpsa_method() samples=10 seconds=60
        
        println("\nBenchmarking particle tracking method...")
        particle_benchmark = @benchmark benchmark_particle_method() samples=10 seconds=60
        
        # Display results
        println("\n" * "="^60)
        println("BENCHMARK RESULTS")
        println("="^60)
        
        println("\nTPSA Method:")
        println("  Minimum time: $(minimum(tpsa_benchmark.times) / 1e6) ms")
        println("  Median time:  $(median(tpsa_benchmark.times) / 1e6) ms")
        println("  Mean time:    $(mean(tpsa_benchmark.times) / 1e6) ms")
        println("  Maximum time: $(maximum(tpsa_benchmark.times) / 1e6) ms")
        println("  Memory allocated: $(tpsa_benchmark.memory) bytes")
        println("  Allocations: $(tpsa_benchmark.allocs)")
        
        println("\nParticle Tracking Method:")
        println("  Minimum time: $(minimum(particle_benchmark.times) / 1e6) ms")
        println("  Median time:  $(median(particle_benchmark.times) / 1e6) ms")
        println("  Mean time:    $(mean(particle_benchmark.times) / 1e6) ms")
        println("  Maximum time: $(maximum(particle_benchmark.times) / 1e6) ms")
        println("  Memory allocated: $(particle_benchmark.memory) bytes")
        println("  Allocations: $(particle_benchmark.allocs)")
        
        # Calculate speedup
        tpsa_median = median(tpsa_benchmark.times)
        particle_median = median(particle_benchmark.times)
        speedup = particle_median / tpsa_median
        
        println("\n" * "="^60)
        println("PERFORMANCE COMPARISON")
        println("="^60)
        println("Speedup (median): $(round(speedup, digits=2))x")
        if speedup > 1
            println("TPSA method is $(round(speedup, digits=2))x FASTER than particle tracking")
        else
            println("Particle tracking is $(round(1/speedup, digits=2))x FASTER than TPSA method")
        end
        
        memory_ratio = particle_benchmark.memory / tpsa_benchmark.memory
        println("Memory usage ratio: $(round(memory_ratio, digits=2))x")
        if memory_ratio > 1
            println("Particle tracking uses $(round(memory_ratio, digits=2))x MORE memory than TPSA")
        else
            println("TPSA uses $(round(1/memory_ratio, digits=2))x MORE memory than particle tracking")
        end
        
        alloc_ratio = particle_benchmark.allocs / tpsa_benchmark.allocs
        println("Allocation ratio: $(round(alloc_ratio, digits=2))x")
        if alloc_ratio > 1
            println("Particle tracking has $(round(alloc_ratio, digits=2))x MORE allocations than TPSA")
        else
            println("TPSA has $(round(1/alloc_ratio, digits=2))x MORE allocations than particle tracking")
        end
        
        println("\n" * "="^60)
        
        return tpsa_benchmark, particle_benchmark
    end
end


begin
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
            
            # Transform FLAME matrix to JuTrack matrix (same as original method)
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
            
            # println("Loaded FLAME data for charge $charge:")
            # println("  Original centroid[6]: $(centroid[6])")
            # println("  Transformed centroid[6]: $(centroid_transformed[6])")
            # println("  Momentum deviation: $momentum_deviation")
        end
        
        return flame_matrices, flame_centroids, settings
    end

    function create_tpsa_beam(centroid)
        tpsa_coords = Vector{CTPS{Float64, 6, 2}}(undef, 6)
        for i in 1:6
            tpsa_coords[i] = CTPS(centroid[i], i, 6,2 )
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

    function propagate_flame_data_TPSA(lattice, flame_matrices, flame_centroids, charges)
        n_elements = length(lattice)
        
        # Initialize storage
        transfer_matrices = Dict{Int, Array{Float64, 3}}()
        beam_matrices = Dict{Int, Array{Float64, 3}}()
        centroids_tpsa = Dict{Int, Matrix{Float64}}()
        
        for charge in charges
            transfer_matrices[charge] = zeros(6, 6, n_elements)
            beam_matrices[charge] = zeros(6, 6, n_elements)
            centroids_tpsa[charge] = zeros(n_elements, 6)
        end
        
        floor_distance = zeros(n_elements)
        
        # Process each charge state
        for charge in charges
            initial_centroid = flame_centroids[charge]
            initial_matrix = flame_matrices[charge]
            
            # println("\nProcessing charge $charge:")
            # println("  Initial centroid: $initial_centroid")
            # println("  Initial matrix determinant: $(det(initial_matrix))")
            
            # Create TPSA beam from FLAME data
            tpsa_coords = create_tpsa_beam(initial_centroid)
            
            # Track through lattice
            for i in eachindex(lattice)
                # Track through element
                pass_TPSA!(lattice[i], tpsa_coords)
                
                # Extract transfer matrix and centroid
                R = extract_transfer_matrix(tpsa_coords)
                transfer_matrices[charge][:, :, i] = R
                centroids_tpsa[charge][i, :] = extract_centroid(tpsa_coords)
                
                # Apply transfer matrix to initial covariance matrix: Σ = R * Σ₀ * R^T
                beam_matrices[charge][:, :, i] = R * initial_matrix * R'
                
                # Calculate floor distance
                if i == 1
                    floor_distance[i] = lattice[i].len
                else
                    floor_distance[i] = floor_distance[i-1] + lattice[i].len
                end
            end
            
            # println("  Final centroid: $(centroids_tpsa[charge][end, :])")
            # println("  Final transfer matrix determinant: $(det(transfer_matrices[charge][:, :, end]))")
        end
        
        return transfer_matrices, beam_matrices, centroids_tpsa, floor_distance
    end

    function extract_beam_properties_from_matrices(beam_matrices, centroids_tpsa, floor_distance, charges)
        
        beam_rms = Dict{Int, Matrix{Float64}}()
        beam_mean = Dict{Int, Matrix{Float64}}()
        beam_centered_rms = Dict{Int, Matrix{Float64}}()
        twi = Dict{Int, Matrix{Float64}}()
        
        n_elements = length(floor_distance)
        
        for charge in charges
            beam_rms[charge] = zeros(n_elements, 6)
            beam_mean[charge] = zeros(n_elements, 6)
            beam_centered_rms[charge] = zeros(n_elements, 6)
            twi[charge] = zeros(n_elements, 9)
            
            for i in 1:n_elements
                # Get covariance matrix at this element
                sigma_matrix = beam_matrices[charge][:, :, i]
                
                # Mean values (centroids)
                beam_mean[charge][i, :] = centroids_tpsa[charge][i, :]
                
                # RMS values are square roots of diagonal elements
                for dim in 1:6
                    beam_rms[charge][i, dim] = sqrt(abs(sigma_matrix[dim, dim]))
                    beam_centered_rms[charge][i, dim] = sqrt(abs(sigma_matrix[dim, dim]))
                end
                
                # Calculate Twiss parameters
                if sigma_matrix[1,1] > 0 && sigma_matrix[2,2] > 0
                    det_x = sigma_matrix[1,1] * sigma_matrix[2,2] - sigma_matrix[1,2]^2
                    if det_x > 0
                        emit_x = sqrt(det_x)
                        beta_x = sigma_matrix[1,1] / emit_x
                        alpha_x = -sigma_matrix[1,2] / emit_x
                        gamma_x = sigma_matrix[2,2] / emit_x
                        
                        twi[charge][i, 1] = beta_x
                        twi[charge][i, 2] = alpha_x
                        twi[charge][i, 3] = emit_x
                    end
                end
                
                if sigma_matrix[3,3] > 0 && sigma_matrix[4,4] > 0
                    det_y = sigma_matrix[3,3] * sigma_matrix[4,4] - sigma_matrix[3,4]^2
                    if det_y > 0
                        emit_y = sqrt(det_y)
                        beta_y = sigma_matrix[3,3] / emit_y
                        alpha_y = -sigma_matrix[3,4] / emit_y
                        gamma_y = sigma_matrix[4,4] / emit_y
                        
                        twi[charge][i, 4] = beta_y
                        twi[charge][i, 5] = alpha_y
                        twi[charge][i, 6] = emit_y
                    end
                end
                
                if sigma_matrix[5,5] > 0 && sigma_matrix[6,6] > 0
                    det_z = sigma_matrix[5,5] * sigma_matrix[6,6] - sigma_matrix[5,6]^2
                    if det_z > 0
                        emit_z = sqrt(det_z)
                        beta_z = sigma_matrix[5,5] / emit_z
                        alpha_z = -sigma_matrix[5,6] / emit_z
                        gamma_z = sigma_matrix[6,6] / emit_z
                        
                        twi[charge][i, 7] = beta_z
                        twi[charge][i, 8] = alpha_z
                        twi[charge][i, 9] = emit_z
                    end
                end
            end
        end
        
        return beam_rms, beam_mean, beam_centered_rms, twi
    end

    function run_simulation_TPSA()
        
        # Load raw FLAME data from YAML
        flame_matrices, flame_centroids, settings = load_flame_data(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        # Load lattice (but not the pre-created beams)
        lattice, _ = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        charges = sort(collect(keys(flame_matrices)))
        
        # Propagate FLAME data using TPSA
        transfer_matrices, beam_matrices, centroids_tpsa, floor_distance = 
            propagate_flame_data_TPSA(lattice, flame_matrices, flame_centroids, charges)
        
        # Extract beam properties from matrices
        beam_rms, beam_mean, beam_centered_rms, twi = 
            extract_beam_properties_from_matrices(beam_matrices, centroids_tpsa, floor_distance, charges)
        
        return lattice, flame_matrices, beam_rms, beam_mean, beam_centered_rms, twi, floor_distance, transfer_matrices
    end

    function compare_with_original_particle_method()
        
        println("\n=== Running Original Particle Tracking Method ===")
        
        # Run original method (creates particle distributions from FLAME data)
        lattice_orig, multibeam_orig = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        # Propagate particle distributions
        function propagate_multibeam_orig(lattice, multibeam::MultiChargeBeam)
            n_elements = length(lattice)
            charges = sort(collect(keys(multibeam.beams)))
            
            beam_rms = Dict{Int, Matrix{Float64}}()
            beam_mean = Dict{Int, Matrix{Float64}}()
            beam_centered_rms = Dict{Int, Matrix{Float64}}()
            twi = Dict{Int, Matrix{Float64}}()
            
            for charge in charges
                beam_rms[charge] = zeros(n_elements, 6)
                beam_mean[charge] = zeros(n_elements, 6)
                beam_centered_rms[charge] = zeros(n_elements, 6)
                twi[charge] = zeros(n_elements, 9)
            end
            
            floor_distance = zeros(n_elements)
            
            for i in eachindex(lattice)
                pass!(lattice[i], multibeam)
                
                for charge in charges
                    beam = multibeam.beams[charge]
                    
                    twi[charge][i, :] .= twiss_beam(beam)
                    means = mean(beam.r, dims=1)
                    beam_mean[charge][i, :] .= vec(means)
                    
                    for dim in 1:6
                        beam_rms[charge][i, dim] = sqrt(mean(beam.r[:,dim].^2))
                        beam_centered_rms[charge][i, dim] = sqrt(mean((beam.r[:,dim] .- means[dim]).^2))
                    end
                end
                
                if i == 1
                    floor_distance[i] = lattice[i].len
                else
                    floor_distance[i] = floor_distance[i-1] + lattice[i].len
                end
            end
            
            return beam_rms, beam_mean, beam_centered_rms, twi, floor_distance
        end
        
        beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, floor_distance_orig = 
            propagate_multibeam_orig(lattice_orig, multibeam_orig)
        
        println("Original particle tracking completed!")
        
        return beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, floor_distance_orig
    end

    function compare_results(beam_rms_tpsa, beam_mean_tpsa, beam_centered_rms_tpsa, twi_tpsa,
                                beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, charges)
        
        println("\n=== Comparing TPSA vs Original Particle Tracking ===")
        
        # Compare RMS values
        println("\nComparing RMS values:")
        for charge in charges
            abs_diff = abs.(beam_centered_rms_tpsa[charge] - beam_centered_rms_orig[charge])
            max_abs_diff = maximum(abs_diff)
            
            rel_diff = zeros(size(beam_centered_rms_orig[charge]))
            for i in eachindex(beam_centered_rms_orig[charge])
                if abs(beam_centered_rms_orig[charge][i]) > 1e-10
                    rel_diff[i] = abs_diff[i] / abs(beam_centered_rms_orig[charge][i])
                end
            end
            
            max_rel_diff = maximum(rel_diff)
            mean_rel_diff = mean(rel_diff[rel_diff .> 0])
            
            println("  Charge $charge:")
            println("    Max absolute difference: $(max_abs_diff)")
            println("    Max relative difference: $(max_rel_diff)")
            println("    Mean relative difference: $(mean_rel_diff)")
        end
        
        # Compare centroids
        println("\nComparing beam centroids:")
        for charge in charges
            abs_diff = abs.(beam_mean_tpsa[charge] - beam_mean_orig[charge])
            max_abs_diff = maximum(abs_diff)
            
            println("  Charge $charge:")
            println("    Max absolute centroid difference: $(max_abs_diff)")
            
            if max_abs_diff > 1e-4
                println("    WARNING: Large centroid differences detected!")
            end
        end
        
        # Compare Twiss parameters
        println("\nComparing Twiss parameters:")
        for charge in charges
            abs_diff = abs.(twi_tpsa[charge] - twi_orig[charge])
            max_abs_diff = maximum(abs_diff)
            
            rel_diff = zeros(size(twi_orig[charge]))
            for i in eachindex(twi_orig[charge])
                if abs(twi_orig[charge][i]) > 1e-10
                    rel_diff[i] = abs_diff[i] / abs(twi_orig[charge][i])
                end
            end
            
            max_rel_diff = maximum(rel_diff[isfinite.(rel_diff)])
            mean_rel_diff = mean(rel_diff[isfinite.(rel_diff) .& (rel_diff .> 0)])
            
            println("  Charge $charge:")
            println("    Max absolute difference: $(max_abs_diff)")
            println("    Max relative difference: $(max_rel_diff)")
            println("    Mean relative difference: $(mean_rel_diff)")
        end
    end
end

# Run the TPSA simulation starting from raw FLAME data
begin
    println("=== Starting TPSA Simulation from FLAME Data ===")
    lattice_tpsa, flame_matrices, beam_rms_tpsa, beam_mean_tpsa, beam_centered_rms_tpsa, twi_tpsa, floor_distance_tpsa, transfer_matrices = run_simulation_TPSA(); # TO precompile
    tpsa_time = @elapsed begin
        lattice_tpsa, flame_matrices, beam_rms_tpsa, beam_mean_tpsa, beam_centered_rms_tpsa, twi_tpsa, floor_distance_tpsa, transfer_matrices = run_simulation_TPSA();
    end
end

begin
    charges = sort(collect(keys(flame_matrices)))
    println("\nTPSA simulation completed in $(round(tpsa_time, digits=3)) seconds")
    println("Number of elements: $(length(lattice_tpsa))")
    println("Charges processed: $charges")
end

# Print final transfer matrix properties
begin
    println("\n=== Final Transfer Matrix Properties ===")
    final_element = length(lattice_tpsa)
    for charge in charges
        final_matrix = transfer_matrices[charge][:, :, final_element]
        det_val = det(final_matrix)
        println("Charge $charge: Final transfer matrix determinant = $det_val")
    end
end

# Symplecticity check
begin
    println("\n=== Symplecticity Check ===")
    S = [0 1 0 0 0 0;
        -1 0 0 0 0 0; 
        0 0 0 1 0 0;
        0 0 -1 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 -1 0]

    for charge in charges
        final_R = transfer_matrices[charge][:, :, end]
        symplectic_error = maximum(abs.(final_R' * S * final_R - S))
        println("Charge $charge: Symplectic error = $symplectic_error")
    end
end

# Generate plots with TPSA suffix
begin
    println("\n=== Generating TPSA Plots ===")
    plts_tpsa = plot_multibeam_data(floor_distance_tpsa, beam_centered_rms_tpsa, twi_tpsa, lattice_tpsa; combined = false);
    display(plts_tpsa["rms"])
    
    # Create a dummy multibeam for plotting 
    multibeam_dummy = MultiChargeBeam(
        Dict(charge => Beam(1e9, 1, 1, charge=Float64(charge)) for charge in charges),
        Dict(charge => 1.0 for charge in charges),
        charges[1]
    )
    
    orbits_tpsa = plot_multicharge_orbits(floor_distance_tpsa, multibeam_dummy, beam_mean_tpsa, save_path="src/demo/FRIB/orbit_plot_TPSA.png")
end

# Export CSV data with TPSA suffix
begin
    println("\n=== Exporting TPSA Data ===")
    export_beam_data_csv(floor_distance_tpsa, beam_mean_tpsa, beam_centered_rms_tpsa, 
                        charges=charges, charge_order=[50, 49, 51], output_dir="src/demo/FRIB/data_TPSA")
end


# # Timing Comparison with original particle tracking
# particle_time = @elapsed begin
#     beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, floor_distance_orig = compare_with_original_particle_method();
# end
# println("Particle tracking completed in $(round(particle_time, digits=3)) seconds")
# println("Speedup: $(round(particle_time/ $(round(tpsa_time, digits=3)), digits=2))x")

#compare TPSA with original
beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, floor_distance_orig = compare_with_original_particle_method();
compare_results(beam_rms_tpsa, beam_mean_tpsa, beam_centered_rms_tpsa, twi_tpsa, beam_rms_orig, beam_mean_orig, beam_centered_rms_orig, twi_orig, charges)

println("\n" * "="^60)
println("DETAILED PERFORMANCE BENCHMARK")
println("="^60)
tpsa_bench, particle_bench = run_performance_comparison()

