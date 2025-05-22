begin
    include("../../JuTrack.jl")
    using .JuTrack
    # using LinearAlgebra
    # using Random
    # using Distributions
    # using CairoMakie
    using Statistics
    # using ProgressMeter
    # using CSV, DataFrames
end

begin
    function propagate_multibeam(lattice, multibeam::MultiChargeBeam)
        n_elements = length(lattice)
        charges = sort(collect(keys(multibeam.beams)))
        
        # Initialize storage for each beam
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
        
        # Propagate through lattice
        for i in eachindex(lattice)
            # Pass all beams through element
            pass!(lattice[i], multibeam)
            
            # Calculate properties for each beam
            for charge in charges
                beam = multibeam.beams[charge]
                
                twi[charge][i, :] .= twiss_beam(beam)
                means = mean(beam.r, dims=1)
                beam_mean[charge][i, :] .= vec(means)
                
                # Store RMS values
                for dim in 1:6
                    beam_rms[charge][i, dim] = sqrt(mean(beam.r[:,dim].^2))
                    beam_centered_rms[charge][i, dim] = sqrt(mean((beam.r[:,dim] .- means[dim]).^2))
                end
            end
            
            # Calculate floor distance
            if i == 1
                floor_distance[i] = lattice[i].len
            else
                floor_distance[i] = floor_distance[i-1] + lattice[i].len
            end
        end
        
        return beam_rms, beam_mean, beam_centered_rms, twi, floor_distance, lattice
    end

    function run_simulation()
        # Load single lattice with all beams
        lattice, multibeam = load_multicharge_lattice(joinpath(pwd(),"src/demo/FRIB/FRIB.yaml"))
        
        # Propagate all beams together
        beam_rms, beam_mean, beam_centered_rms, twi, floor_distance, lattice = 
            propagate_multibeam(lattice, multibeam)
        
        return lattice, multibeam, beam_rms, beam_mean, beam_centered_rms, twi, floor_distance
    end
end

lattice, multibeam, beam_rms, beam_mean, beam_centered_rms, twi, floor_distance = run_simulation();


# Plot results 
begin
plts = plot_multibeam_data(floor_distance, beam_centered_rms, twi, lattice; combined = false);
display(plts["rms"])
orbits = plot_multicharge_orbits(floor_distance, multibeam, beam_mean, save_path="src/demo/FRIB/orbit_plot.png")
end

# Access individual beams if needed
beam50 = multibeam.beams[50]
beam49 = multibeam.beams[49]
beam51 = multibeam.beams[51]


export_beam_data_csv(floor_distance, beam_mean, beam_centered_rms, 
                    charge_order=[50, 49, 51])

function print_lattice(lattice)
    println("Lattice Elements:")
    for i in eachindex(lattice)
        println("Element $(i): $(typeof(lattice[i])) - $(lattice[i].name)")
    end
end

print_lattice(lattice)



function get_flame_center(beam_mean)
    scaling =  [0.001, 1.0, 0.001, 1.0, 0.35242784364090085, 0.0024413881418195262]
    new_center_for_flame = beam_mean[1,:] ./ scaling
    push!(new_center_for_flame, 1.0)
    println(new_center_for_flame)
end

# get_flame_center(beam49_mean)
# get_flame_center(beam51_mean)