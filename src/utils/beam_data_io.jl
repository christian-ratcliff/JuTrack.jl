using CSV, DataFrames

# Universal CSV export function - works for any charge states
function export_beam_data_csv(floor_distance, beam_mean_data=nothing, beam_rms_data=nothing; 
                             charges=nothing, charge_order=nothing, output_dir="src/demo/FRIB/data")
    """
    Export beam data (mean and RMS) to CSV files.
    
    Args:
        floor_distance: Array of s positions
        beam_mean_data: Either Dict{Int, Matrix} (MultiChargeBeam) or nothing (global vars)
        beam_rms_data: Either Dict{Int, Matrix} (MultiChargeBeam) or nothing (global vars)  
        charges: List of charge states to export (auto-detected if not provided)
        output_dir: Directory to save CSV files
    """
    
    # Auto-detect charges if not provided
    if charges === nothing
        if isa(beam_mean_data, Dict)
            charges = sort(collect(keys(beam_mean_data)))
        elseif isa(beam_rms_data, Dict)
            charges = sort(collect(keys(beam_rms_data)))
        else
            # Try to find global variables
            charges = []
            for i in 1:200  # Check reasonable range
                try
                    eval(Symbol("beam$(i)_mean"))
                    push!(charges, i)
                catch
                    # Variable doesn't exist, skip
                end
            end
            if isempty(charges)
                error("No beam data found. Either provide data dictionaries or ensure beam variables exist.")
            end
        end
    end
    
    println("Exporting data for charges: $charges")
    
    # Create output directory if it doesn't exist
    mkpath(output_dir)

    reference_charge = charge_order[1]  # Use the first charge in the provided order as reference
    
    # Create mapping with reference beam as beam1
    charge_to_beam = Dict{Int, Int}()
    charge_to_beam[reference_charge] = 1
    
    # Add other charges as beam2, beam3, etc.
    beam_counter = 2
    for charge in sort(charges)
        if charge != reference_charge
            charge_to_beam[charge] = beam_counter
            beam_counter += 1
        end
    end
    
    println("Reference charge: $reference_charge (beam1)")
    
    # Initialize DataFrames
    df_centered_rms = DataFrame(s = floor_distance)
    df_mean = DataFrame(s = floor_distance)
    
    # Process each charge state in the specified order
    for charge in charge_order
        beam_num = charge_to_beam[charge]
        
        # Get data based on input type
        if isa(beam_mean_data, Dict) && isa(beam_rms_data, Dict)
            # MultiChargeBeam case
            if haskey(beam_mean_data, charge) && haskey(beam_rms_data, charge)
                mean_data = beam_mean_data[charge]
                rms_data = beam_rms_data[charge]
                
                # Check and fix dimensions if needed
                if size(mean_data, 1) != length(floor_distance)
                    if size(mean_data, 2) == length(floor_distance)
                        mean_data = mean_data'
                    else
                        @warn "Cannot fix dimension mismatch for charge $charge mean data, skipping"
                        continue
                    end
                end
                
                if size(rms_data, 1) != length(floor_distance)
                    if size(rms_data, 2) == length(floor_distance)
                        rms_data = rms_data'
                    else
                        @warn "Cannot fix dimension mismatch for charge $charge RMS data, skipping"
                        continue
                    end
                end
            else
                @warn "Charge $charge not found in beam data dictionaries, skipping"
                continue
            end
        else
            # Legacy case: use global variables
            try
                mean_data = eval(Symbol("beam$(charge)_mean"))
                rms_data = eval(Symbol("beam$(charge)_centered_rms"))
            catch e
                @warn "Could not find variables for charge $charge: $e"
                continue
            end
        end
        
        # Add RMS data columns
        df_centered_rms[!, "beam$(beam_num)_charge$(charge)_x_rms"] = rms_data[:,1]
        df_centered_rms[!, "beam$(beam_num)_charge$(charge)_xp_rms"] = rms_data[:,2] 
        df_centered_rms[!, "beam$(beam_num)_charge$(charge)_y_rms"] = rms_data[:,3]
        df_centered_rms[!, "beam$(beam_num)_charge$(charge)_yp_rms"] = rms_data[:,4]
        
        # Add mean data columns
        df_mean[!, "beam$(beam_num)_charge$(charge)_x_mean"] = mean_data[:,1]
        df_mean[!, "beam$(beam_num)_charge$(charge)_xp_mean"] = mean_data[:,2]
        df_mean[!, "beam$(beam_num)_charge$(charge)_y_mean"] = mean_data[:,3]
        df_mean[!, "beam$(beam_num)_charge$(charge)_yp_mean"] = mean_data[:,4]
    end
    
    # Write CSV files
    rms_path = joinpath(output_dir, "all_beams_centered_rms.csv")
    mean_path = joinpath(output_dir, "all_beams_mean.csv")
    
    CSV.write(rms_path, df_centered_rms)
    CSV.write(mean_path, df_mean)
    
    println("CSV files saved:")
    println("  RMS data: $rms_path")
    println("  Mean data: $mean_path")
    println("  Charge to beam mapping: $charge_to_beam")
    
    return nothing
end