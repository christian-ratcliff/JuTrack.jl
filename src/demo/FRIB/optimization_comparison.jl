begin
    include("../../JuTrack.jl")
    using .JuTrack
    using LinearAlgebra
    using Random
    using Statistics
    using Printf
    using Enzyme
end





# Load the lattice and beam data
function setup_optimization_problem(error_factor=0.5)
    # Load lattice and multibeam
    lattice, multibeam = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
    
    # Apply error to initial beam conditions
    Random.seed!(42)
    perturbed_multibeam = perturb_beam_conditions(multibeam, error_factor)
    
    # Identify quadrupoles for optimization
    quad_indices = identify_optimization_quads(lattice)
    
    return lattice, perturbed_multibeam, quad_indices
end


function setup_optimization_elements(lattice, quad_indices)
    element_info = []
    for idx in quad_indices
        elem = lattice[idx]
        push!(element_info, (name=elem.name, len=elem.len))
    end
    return element_info
end


# Perturb the initial beam conditions
function perturb_beam_conditions(multibeam::MultiChargeBeam, error_factor::Float64)
    perturbed_beams = Dict{Int, Beam}()
    
    for (charge, beam) in multibeam.beams
        # Create a copy of the beam
        new_beam = Beam(beam)
        
        # Get initial 2nd moment matrix
        moment2nd = copy(beam.moment2nd)
        
        # Apply random error to each matrix element
        for i in 1:6
            for j in i:6
                error = 1.0 + error_factor * (2*rand() - 1)  # Random error between ±error_factor
                moment2nd[i,j] *= error
                moment2nd[j,i] = moment2nd[i,j]  # Keep symmetric
            end
        end
        
        # Ensure the matrix is positive definite
        eigenvals, eigenvecs = eigen(moment2nd)
        eigenvals = max.(eigenvals, 1e-10)  # Ensure positive eigenvalues
        moment2nd = eigenvecs * Diagonal(eigenvals) * eigenvecs'
        
        # Generate new particle distribution from perturbed matrix
        nparticles = beam.np
        dis = randn(nparticles, 6)
        
        # Transform to match perturbed covariance
        lam, u = eigen(moment2nd)
        transformation = u * Diagonal(sqrt.(lam)) * u'
        dis = dis * transformation'
        
        # Add centroid
        for i in 1:6
            dis[:, i] .+= beam.centroid[i]
        end
        
        new_beam.r = dis
        new_beam.moment2nd = moment2nd
        
        # Recalculate beam properties
        get_emittance!(new_beam)
        get_centroid!(new_beam)
        
        perturbed_beams[charge] = new_beam
    end
    
    return MultiChargeBeam(perturbed_beams, multibeam.brho_values, multibeam.reference_charge)
end

# Identify quadrupoles to use for optimization
function identify_optimization_quads(lattice)
    quad_indices = Int[]
    bend_start = 0
    bend_end = 0
    
    # Find bending section boundaries
    for (i, elem) in enumerate(lattice)
        if elem isa SBEND && bend_start == 0
            bend_start = i
        elseif elem isa SBEND && bend_start > 0
            bend_end = i
        end
    end
    
    # Select quadrupoles around bending section
    # Before bending: look for quads in range [bend_start-200, bend_start-20]
    # After bending: look for quads in range [bend_end+20, bend_end+200]
    
    for (i, elem) in enumerate(lattice)
        if elem isa KQUAD
            # Before bending section
            if (bend_start - 200) <= i <= (bend_start - 20)
                push!(quad_indices, i)
            # After bending section  
            elseif (bend_end + 20) <= i <= (bend_end + 200)
                push!(quad_indices, i)
            end
        end
    end
    
    # Limit to 10-15 quadrupoles
    # if length(quad_indices) > 10
    #     # Select evenly spaced quads
    #     step = length(quad_indices) ÷ 10
    #     quad_indices = quad_indices[1:step:end][1:10]
    # end
    
    println("Found $(length(quad_indices)) possible quadrupoles for optimization")
    println("Quadrupole indices: ", quad_indices)
    
    return quad_indices
end

# Define locations for beam size evaluation
function get_evaluation_points(lattice)
    eval_points = Dict{String, Vector{Int}}()
    
    # Find first straight section (before first bend)
    first_straight_end = 0
    for (i, elem) in enumerate(lattice)
        if elem isa SBEND
            first_straight_end = i - 1
            break
        end
    end
    eval_points["first_straight"] = collect(1:10:first_straight_end)
    
    # Find after bend section
    last_bend = 0
    for (i, elem) in enumerate(lattice)
        if elem isa SBEND
            last_bend = i
        end
    end
    eval_points["after_bend"] = collect(last_bend+1:10:length(lattice)-1)
    
    # End point
    eval_points["end"] = [length(lattice)]
    
    return eval_points
end

# Perturb FLAME data for TPSA
function perturb_flame_data(flame_matrices, flame_centroids, error_factor::Float64)
    Random.seed!(42) 
    
    perturbed_matrices = Dict{Int, Matrix{Float64}}()
    perturbed_centroids = Dict{Int, Vector{Float64}}()
    
    for (charge, matrix) in flame_matrices
        perturbed_matrix = copy(matrix)
        for i in 1:6
            for j in i:6
                error = 1.0 + error_factor * (2*rand() - 1)
                perturbed_matrix[i,j] *= error
                perturbed_matrix[j,i] = perturbed_matrix[i,j]  # Keep symmetric
            end
        end
        
        # Ensure positive definite
        eigenvals, eigenvecs = eigen(perturbed_matrix)
        eigenvals = max.(eigenvals, 1e-10)
        perturbed_matrix = eigenvecs * Diagonal(eigenvals) * eigenvecs'
        
        # Perturb centroid slightly
        perturbed_centroid = copy(flame_centroids[charge])
        for i in 1:6
            centroid_error = error_factor * 0.1 * (2*rand() - 1)  
            perturbed_centroid[i] *= (1.0 + centroid_error)
        end
        
        perturbed_matrices[charge] = perturbed_matrix
        perturbed_centroids[charge] = perturbed_centroid
    end
    
    return perturbed_matrices, perturbed_centroids
end

# TPSA cost function using segments
# function tpsa_cost(k_values::Vector{Float64}, lattice, flame_matrices, flame_centroids, 
#                    quad_indices::Vector{Int}, element_info::Vector, 
#                    eval_points::Dict, charges::Vector{Int})
    
#     # Create changed elements
#     changed_ele = []
#     for (i, info) in enumerate(element_info)
#         new_quad = KQUAD(name=info.name, len=info.len, k1=k_values[i])
#         push!(changed_ele, new_quad)
#     end
    
#     cost = 0.0
    
#     # Find segment boundaries
#     first_bend_start = findfirst(elem -> elem isa SBEND, lattice)
#     last_bend_end = findlast(elem -> elem isa SBEND, lattice)
    
#     # Process each charge state
#     for charge in charges
#         # Create TPSA beam
#         tpsa_coords = create_tpsa_beam(flame_centroids[charge])
#         initial_sigma = flame_matrices[charge]
        
#         # SEGMENT 1: Track from start to end of first section (before bend)
#         segment1 = lattice[1:first_bend_start-1]
#         segment1_quad_indices = Int[]
#         segment1_changed_ele = []
        
#         for (i, global_idx) in enumerate(quad_indices)
#             if 1 <= global_idx <= first_bend_start-1
#                 push!(segment1_quad_indices, global_idx)
#                 push!(segment1_changed_ele, changed_ele[i])
#             end
#         end
        
#         ADlinepass_TPSA!(segment1, tpsa_coords, segment1_quad_indices, segment1_changed_ele)
        
#         # CHECK 1: End of first section - beam size should be < 3mm
#         R1 = extract_transfer_matrix(tpsa_coords)
#         sigma1 = R1 * initial_sigma * R1'
#         x_rms1 = sqrt(abs(sigma1[1, 1])) * 1000
#         y_rms1 = sqrt(abs(sigma1[3, 3])) * 1000
#         # beam_size_first_section = max(x_rms1, y_rms1)
        
#         # UPDATE COST 1: First straight section constraint
#         # if beam_size_first_section > 3.0
#         #     cost += (beam_size_first_section - 3.0)^2 * 100
#         # end

#         cost += ((x_rms1 - 3.)^10 + (y_rms1 - 3.)^10) * 100
        
#         # SEGMENT 2: Track through bend region
#         segment2 = lattice[first_bend_start:last_bend_end]
#         segment2_quad_indices = Int[]
#         segment2_changed_ele = []
        
#         for (i, global_idx) in enumerate(quad_indices)
#             if first_bend_start <= global_idx <= last_bend_end
#                 local_idx = global_idx - first_bend_start + 1
#                 push!(segment2_quad_indices, local_idx)
#                 push!(segment2_changed_ele, changed_ele[i])
#             end
#         end
        
#         ADlinepass_TPSA!(segment2, tpsa_coords, segment2_quad_indices, segment2_changed_ele)
        
#         # CHECK 2: After bend - beam size should be < 3.5mm
#         R2 = extract_transfer_matrix(tpsa_coords)
#         sigma2 = R2 * initial_sigma * R2'
#         x_rms2 = sqrt(abs(sigma2[1, 1])) * 1000
#         y_rms2 = sqrt(abs(sigma2[3, 3])) * 1000
#         # beam_size_after_bend = max(x_rms2, y_rms2)
        
#         # UPDATE COST 2: After bend constraint
#         # if beam_size_after_bend > 3.5
#         #     cost += (beam_size_after_bend - 3.5)^2 * 100
#         # end

#         cost += ((x_rms2 - 3.5)^10 + (y_rms2 - 3.5)^10) * 100
        
#         # SEGMENT 3: Track from after bend to very end
#         segment3 = lattice[last_bend_end+1:end]
#         segment3_quad_indices = Int[]
#         segment3_changed_ele = []
        
#         for (i, global_idx) in enumerate(quad_indices)
#             if global_idx > last_bend_end
#                 local_idx = global_idx - last_bend_end
#                 push!(segment3_quad_indices, local_idx)
#                 push!(segment3_changed_ele, changed_ele[i])
#             end
#         end
        
#         ADlinepass_TPSA!(segment3, tpsa_coords, segment3_quad_indices, segment3_changed_ele)
        
#         # CHECK 3: Very end - beam size should be ~0.3mm
#         R3 = extract_transfer_matrix(tpsa_coords)
#         sigma3 = R3 * initial_sigma * R3'
#         x_rms3 = sqrt(abs(sigma3[1, 1])) * 1000
#         y_rms3 = sqrt(abs(sigma3[3, 3])) * 1000
#         # beam_size_final = max(x_rms3, y_rms3)
        
#         # UPDATE COST 3: Final beam size constraint
#         # cost += (beam_size_final - 0.3)^2 * 1000

#         cost += ((x_rms3 - 0.3)^10 + (y_rms3 - 0.3)^10) * 100
        
#     end
    
#     return cost
# end

function setup_segments(quad_indices, lattice)
    first_bend_start = findfirst(elem -> elem isa SBEND, lattice)
    last_bend_end = findlast(elem -> elem isa SBEND, lattice)
    
    # Pre-compute which quads belong to each segment
    seg1_indices = Int[]
    seg1_positions = Int[]
    seg2_indices = Int[]  
    seg2_positions = Int[]
    seg3_indices = Int[]
    seg3_positions = Int[]
    
    for (i, global_idx) in enumerate(quad_indices)
        if 1 <= global_idx <= first_bend_start-1
            push!(seg1_indices, global_idx)
            push!(seg1_positions, i)
        elseif first_bend_start <= global_idx <= last_bend_end
            push!(seg2_indices, global_idx - first_bend_start + 1)
            push!(seg2_positions, i)
        elseif global_idx > last_bend_end
            push!(seg3_indices, global_idx - last_bend_end)
            push!(seg3_positions, i)
        end
    end
    
    return (first_bend_start, last_bend_end, 
            seg1_indices, seg1_positions,
            seg2_indices, seg2_positions, 
            seg3_indices, seg3_positions)
end

# TPSA cost function without the ifs? Not sure if this is working
function tpsa_cost(k_values::Vector{Float64}, lattice, flame_matrices, flame_centroids, 
                   quad_indices::Vector{Int}, element_info::Vector, 
                   eval_points::Dict, charges::Vector{Int}, segment_info)
    
    # Unpack pre-computed segment info
    first_bend_start, last_bend_end, 
    seg1_indices, seg1_positions,
    seg2_indices, seg2_positions, 
    seg3_indices, seg3_positions = segment_info
    
    # Create changed elements for each segment
    seg1_changed_ele = []
    for pos in seg1_positions
        info = element_info[pos]
        push!(seg1_changed_ele, KQUAD(name=info.name, len=info.len, k1=k_values[pos]))
    end
    
    seg2_changed_ele = []
    for pos in seg2_positions
        info = element_info[pos]
        push!(seg2_changed_ele, KQUAD(name=info.name, len=info.len, k1=k_values[pos]))
    end
    
    seg3_changed_ele = []
    for pos in seg3_positions
        info = element_info[pos]
        push!(seg3_changed_ele, KQUAD(name=info.name, len=info.len, k1=k_values[pos]))
    end
    
    cost = 0.0
    
    # Pre-defined segments
    segment1 = lattice[1:first_bend_start-1]
    segment2 = lattice[first_bend_start:last_bend_end]
    segment3 = lattice[last_bend_end+1:end]
    
    for charge in charges
        # Create TPSA beam
        tpsa_coords = create_tpsa_beam(flame_centroids[charge])
        initial_sigma = flame_matrices[charge]
        
        # SEGMENT 1
        ADlinepass_TPSA!(segment1, tpsa_coords, seg1_indices, seg1_changed_ele)
        
        R1 = extract_transfer_matrix(tpsa_coords)
        sigma1 = R1 * initial_sigma * R1'
        x_rms1 = sqrt(abs(sigma1[1, 1])) 
        y_rms1 = sqrt(abs(sigma1[3, 3])) 
        
        cost += (log(1 + exp(10.0*(x_rms1 - 0.003))) + log(1 + exp(10.0*(y_rms1 - 0.003)))) * 1
        
        # SEGMENT 2
        ADlinepass_TPSA!(segment2, tpsa_coords, seg2_indices, seg2_changed_ele)
        
        R2 = extract_transfer_matrix(tpsa_coords)
        sigma2 = R2 * initial_sigma * R2'
        x_rms2 = sqrt(abs(sigma2[1, 1])) 
        y_rms2 = sqrt(abs(sigma2[3, 3])) 
        
        cost += (log(1 + exp(10.0*(x_rms2 - 0.0035))) + log(1 + exp(10.0*(y_rms2 - 0.0035)))) * 1
        
        # SEGMENT 3
        ADlinepass_TPSA!(segment3, tpsa_coords, seg3_indices, seg3_changed_ele)

        # If only want to look at end, comment out above stuff, and use this
        # ADlinepass_TPSA!(lattice, tpsa_coords, [seg1_indices; seg2_indices; seg3_indices], [seg1_changed_ele; seg2_changed_ele; seg3_changed_ele])
        
        R3 = extract_transfer_matrix(tpsa_coords)
        sigma3 = R3 * initial_sigma * R3'
        x_rms3 = sqrt(abs(sigma3[1, 1])) 
        y_rms3 = sqrt(abs(sigma3[3, 3])) 
        
        cost += ((x_rms3 - 0.0003)^2 + (y_rms3 - 0.0003)^2) * 10
    end
    
    return cost
end

# TPSA optimization function
function optimize_tpsa(lattice, flame_matrices, flame_centroids; 
                      error_factor=0.5, max_iterations=50, step=5e-6
                      , tolerance=1e-6)
    
    println("=== TPSA Optimization with Segments ===")
    
    # Perturb FLAME data
    # perturbed_matrices, perturbed_centroids = perturb_flame_data(flame_matrices, flame_centroids, error_factor)
    perturbed_matrices, perturbed_centroids, _ = load_flame_data(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml")) #After we changed to perturbing k values
    
    # Setup optimization
    quad_indices = identify_optimization_quads(lattice)
    quad_indices = quad_indices[end-9:end]
    element_info = setup_optimization_elements(lattice, quad_indices)
    eval_points = get_evaluation_points(lattice)
    
    charges = sort(collect(keys(perturbed_matrices)))
    k_values = [lattice[i].k1 for i in quad_indices]
    k_values .= k_values .* (1.0 .+ 0.15 .* randn(length(k_values)))
    k_values_original = copy(k_values)

    segment_info = setup_segments(quad_indices, lattice)
    
    # Find and display segment boundaries
    first_bend_start = findfirst(elem -> elem isa SBEND, lattice)
    last_bend_end = findlast(elem -> elem isa SBEND, lattice)
    
    println("Segment boundaries:")
    println("  Segment 1 (first straight): 1 to $(first_bend_start-1)")
    println("  Segment 2 (bend region): $first_bend_start to $last_bend_end") 
    println("  Segment 3 (after bend): $(last_bend_end+1) to $(length(lattice))")
    println("Selected $(length(quad_indices)) quadrupoles: $quad_indices")
    println("Charges: $charges")
    
    # Optimization loop
    cost_history = Float64[]
    
    for iter in 1:max_iterations
        # Calculate current cost


        current_cost = tpsa_cost(k_values, lattice, perturbed_matrices, perturbed_centroids,
                                quad_indices, element_info, eval_points, charges, segment_info)
        push!(cost_history, current_cost)
        
        # Calculate gradients
        gradients = zeros(length(k_values))
        for i in 1:length(k_values)
            direction = zeros(length(k_values))
            direction[i] = 1.0
            grad = autodiff(ForwardWithPrimal, tpsa_cost, Duplicated(k_values, direction),
                           Const(lattice), Const(perturbed_matrices), Const(perturbed_centroids),
                           Const(quad_indices), Const(element_info), Const(eval_points), Const(charges), Const(segment_info))
            gradients[i] = grad[1]
        end
        println("gradients: $gradients")
        # Update k values
        k_values_new = k_values .- step .* gradients
        
        if iter % 1 == 0 || iter == 1
            println("Iteration $iter: Cost = $(round(current_cost, digits=6)), |grad| = $(round(norm(gradients), digits=6)), tolerance = $(round(norm(k_values_new - k_values), digits=4)) ")
            println("k values originally: $k_values_original")
            println("k values now: $k_values_new")
            # println("difference: $(k_values_new - k_values)" )
        end
        
        if norm(k_values_new - k_values) < tolerance
            println("Converged after $iter iterations")
            break
        end
        
        k_values = k_values_new
    end
    
    println("\nTPSA Optimization completed:")
    println("Final cost: $(round(cost_history[end], digits=6))")
    
    return k_values, cost_history, quad_indices, perturbed_matrices, perturbed_centroids
end

# Verification function 
function verify_tpsa_results(lattice, perturbed_matrices, perturbed_centroids, 
                            quad_indices, k_optimal, element_info)
    println("\n=== Verifying TPSA Results ===")
    
    charges = sort(collect(keys(perturbed_matrices)))
    
    # Create optimized elements
    changed_ele = []
    for (i, info) in enumerate(element_info)
        new_quad = KQUAD(name=info.name, len=info.len, k1=k_optimal[i])
        push!(changed_ele, new_quad)
    end
    
    # Find segment boundaries
    first_bend_start = findfirst(elem -> elem isa SBEND, lattice)
    last_bend_end = findlast(elem -> elem isa SBEND, lattice)
    
    for charge in charges
        println("\nCharge $charge:")
        
        # Create TPSA beam
        tpsa_coords = create_tpsa_beam(perturbed_centroids[charge])
        initial_sigma = perturbed_matrices[charge]
        
        # Track through entire lattice with optimized quads
        ADlinepass_TPSA!(lattice, tpsa_coords, quad_indices, changed_ele)
        
        # Calculate final beam size
        R_final = extract_transfer_matrix(tpsa_coords)
        final_sigma = R_final * initial_sigma * R_final'
        
        x_rms_final = sqrt(abs(final_sigma[1, 1])) * 1000
        y_rms_final = sqrt(abs(final_sigma[3, 3])) * 1000
        final_size = max(x_rms_final, y_rms_final)
        
        println("  Final beam size: $(round(final_size, digits=3)) mm (target ~0.3 mm)")
        println("  X RMS: $(round(x_rms_final, digits=3)) mm")
        println("  Y RMS: $(round(y_rms_final, digits=3)) mm")
    end
end

# Main function
function run_tpsa_optimization()
    println("Loading FLAME data and lattice...")
    
    # Load data
    flame_matrices, flame_centroids, _ = load_flame_data(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
    lattice, _ = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
    
    println("Loaded $(length(lattice)) lattice elements")
    
    # Run optimization
    t_start = time()
    k_optimal, cost_history, quad_indices, perturbed_matrices, perturbed_centroids = 
        optimize_tpsa(lattice, flame_matrices, flame_centroids,
                     error_factor=0.5, max_iterations=100, step=5e-3, tolerance=1e-6)
    t_elapsed = time() - t_start
    
    println("\nOptimization completed in $(round(t_elapsed, digits=2)) seconds")
    
    # Verify results
    element_info = setup_optimization_elements(lattice, quad_indices)
    verify_tpsa_results(lattice, perturbed_matrices, perturbed_centroids, 
                       quad_indices, k_optimal, element_info)
    
    return k_optimal, cost_history, quad_indices
end

# Run the TPSA optimization
k_optimal, cost_history, quad_indices = run_tpsa_optimization()


