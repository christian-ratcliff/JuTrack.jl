using CairoMakie

# Universal orbit plotting function - works for any charge states
function plot_orbits(floor_distance, beam_data=nothing; charges=nothing, save_path=nothing)
    """
    Plot orbit data for multiple charge states.
    
    Args:
        floor_distance: Array of s positions
        beam_data: Either Dict{Int, Matrix} (from MultiChargeBeam) or nothing (uses global variables)
        charges: List of charge states to plot (auto-detected if not provided)
        save_path: Optional path to save the plot
    """
    
    # Auto-detect charges if not provided
    if charges === nothing
        if isa(beam_data, Dict)
            charges = sort(collect(keys(beam_data)))
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
                error("No beam data found. Either provide beam_data dictionary or ensure beam variables exist.")
            end
        end
    end
    
    println("Plotting charges: $charges")
    
    
    fig = Figure(size=(2000, 1200))
    ax1 = Axis(fig[1, 1], title="X Position Orbits", xlabel="s [m]", ylabel="[mm]")
    ax2 = Axis(fig[2, 1], title="Y Position Orbits", xlabel="s [m]", ylabel="[mm]")
    ax3 = Axis(fig[1, 2], title="X Momentum Orbits", xlabel="s [m]")
    ax4 = Axis(fig[2, 2], title="Y Momentum Orbits", xlabel="s [m]")
    
    # Generate colors automatically
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    
    legend_elements = []
    legend_labels = []
    
    for (idx, charge) in enumerate(charges)
        color = colors[mod1(idx, length(colors))]  # Cycle through colors
        
        # Get mean data based on input type
        if isa(beam_data, Dict)
            # MultiChargeBeam case: data is in dictionaries
            if haskey(beam_data, charge)
                mean_data = beam_data[charge]
            else
                @warn "Charge $charge not found in beam_data, skipping"
                continue
            end
        else
            # Legacy case: use global variables
            try
                mean_data = eval(Symbol("beam$(charge)_mean"))
            catch
                @warn "Variable beam$(charge)_mean not found, skipping charge $charge"
                continue
            end
        end
        
        # Check dimensions and fix if needed
        if size(mean_data, 1) != length(floor_distance)
            @warn "Dimension mismatch for charge $charge: floor_distance=$(length(floor_distance)), mean_data=$(size(mean_data))"
            if size(mean_data, 2) == length(floor_distance)
                # Data is transposed, fix it
                mean_data = mean_data'
                println("Transposed data for charge $charge, new size: ", size(mean_data))
            else
                @warn "Cannot fix dimension mismatch for charge $charge, skipping"
                continue
            end
        end
        
        x_data = mean_data[:,1] .* 1000
        y_data = mean_data[:,3] .* 1000
        xp_data = mean_data[:,2]
        yp_data = mean_data[:,4]
        
        # Plot position orbits
        xline = lines!(ax1, floor_distance, x_data, 
                    label="x, $charge", color=color)
        lines!(ax2, floor_distance, y_data, 
            label="y, $charge", color=color)
        
        # Plot momentum orbits  
        lines!(ax3, floor_distance, xp_data, 
               label="x', $charge", color=color)
        lines!(ax4, floor_distance, yp_data, 
               label="y', $charge", color=color)
        
        # Store for legend (use x-position line)
        push!(legend_elements, xline)
        push!(legend_labels, string(charge))
    end
    
    # Add legend
    Legend(fig[1,3], legend_elements, legend_labels, "Charge States")
    
    # display(fig)
    
    # Save if path provided
    if save_path !== nothing
        save(save_path, fig)
        println("Orbit plot saved to: $save_path")
    end
    
    return fig
end

# Helper function to get available charges from any beam data source
function get_available_charges(beam_data=nothing)
    """Get list of available charge states from beam data or global variables."""
    if isa(beam_data, Dict)
        return sort(collect(keys(beam_data)))
    elseif isa(beam_data, MultiChargeBeam)
        return sort(collect(keys(beam_data.beams)))
    else
        # Search global variables
        charges = []
        for i in 1:200
            try
                eval(Symbol("beam$(i)_mean"))
                push!(charges, i)
            catch
                # Continue searching
            end
        end
        return charges
    end
end

# Specific function for MultiChargeBeam objects
function plot_multicharge_orbits(floor_distance, multibeam::MultiChargeBeam, beam_mean; save_path=nothing)
    charges = sort(collect(keys(multibeam.beams)))
    return plot_orbits(floor_distance, beam_mean, charges=charges, save_path=save_path)
end


function plot_multibeam_data(floor_length, beam_rms, twi, lattice; 
                            charges=nothing, mass_number=124, element_symbol="Xe", combined=true)
    """
    Generalized plotting function for multiple charge states.
    
    Args:
        floor_length: Array of s positions
        beam_rms: Dict{Int, Matrix} containing RMS data for each charge
        twi: Dict{Int, Matrix} containing Twiss data for each charge  
        lattice: Array of lattice elements
        charges: Array of charge states to plot (auto-detected if not provided)
        mass_number: Mass number for labels (default 124)
        element_symbol: Element symbol for labels (default "Xe")
        combined: If true, create single combined plot; if false, return separate plots
    """
    
    # Auto-detect charges if not provided
    if charges === nothing
        charges = sort(collect(keys(beam_rms)))
    end
    
    println("Plotting charges: $charges")
    
    # Generate colors automatically (cycle through if more charges than colors)
    base_colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
    beam_colors = [base_colors[mod1(i, length(base_colors))] for i in 1:length(charges)]
    
    # Generate labels automatically
    beam_labels = [L"^{%$(mass_number)}%$(element_symbol)^{%$(charge)+}" for charge in charges]
    
    # Extract data for each charge state
    beam_rms_data = []
    twi_data = []
    
    for charge in charges
        if haskey(beam_rms, charge) && haskey(twi, charge)
            push!(beam_rms_data, beam_rms[charge])
            push!(twi_data, twi[charge])
        else
            error("Charge $charge not found in beam_rms or twi data")
        end
    end
    
    if combined != true    
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
        function create_plot(top_data_list, bottom_data_list, top_label, bottom_label, title)
            
            f = Figure(size = (1200, 400))
            gl = f[1, 1] = GridLayout()
            legend_layout = f[1, 2] = GridLayout()
            
            # Top plot
            ax_top = Axis(gl[1, 1], xlabel = "", ylabel = top_label)
            for (i, data) in enumerate(top_data_list)
                lines!(ax_top, floor_length, data, color = beam_colors[i], linewidth = 2)
            end
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
                ele_width = ele.len
                
                # Draw rectangle for element based on element type and properties
                if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
                    # Use the original quad height
                    quad_height = element_heights["KQUAD"]
                    drift_height = element_heights["DRIFT"]
                    
                    # Calculate offset to align with other elements
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
                
                # Update position
                curr_pos += ele_width
            end
            
            hidexdecorations!(floor_ax)
            
            # Bottom plot
            ax_bottom = Axis(gl[3, 1], xlabel = "z [m]", ylabel = bottom_label, yreversed = true)
            for (i, data) in enumerate(bottom_data_list)
                lines!(ax_bottom, floor_length, data, color = beam_colors[i], linewidth = 2)
            end
            
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
            
            return f
        end
        
        # Create the four plots with dynamic data
        rms_plot = create_plot(
            [data[:,1] .* 1e3 for data in beam_rms_data],  # X RMS data for all charges
            [data[:,3] .* 1e3 for data in beam_rms_data],  # Y RMS data for all charges
            "RMS X [mm]",
            "RMS Y [mm]",
            "RMS Plot"
        )
        
        beta_plot = create_plot(
            [data[:,1] for data in twi_data],  # Beta X data for all charges
            [data[:,4] for data in twi_data],  # Beta Y data for all charges
            L"\beta_x",
            L"\beta_y",
            "Beta Functions"
        )
        
        alpha_plot = create_plot(
            [data[:,2] for data in twi_data],  # Alpha X data for all charges
            [data[:,5] for data in twi_data],  # Alpha Y data for all charges
            L"\alpha_x",
            L"\alpha_y",
            "Alpha Functions"
        )
        
        emittance_plot = create_plot(
            [data[:,3] for data in twi_data],  # Emittance X data for all charges
            [data[:,6] for data in twi_data],  # Emittance Y data for all charges
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
        # Create one large combined figure
        f = Figure(size = (1200, 1200))
        
        # Create the main layout
        gl = f[1, 1] = GridLayout()
        legend_layout = f[1, 2] = GridLayout()
        
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
        for (i, data) in enumerate(beam_rms_data)
            lines!(ax_rms_x, floor_length, data[:,1] .* 1e3, color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_rms_x)
        
        # 2. RMS Y
        ax_rms_y = Axis(gl[2, 1], xlabel = "", ylabel = "RMS Y [mm]")
        for (i, data) in enumerate(beam_rms_data)
            lines!(ax_rms_y, floor_length, data[:,3] .* 1e3, color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_rms_y)
        
        # 3. Beta X
        ax_beta_x = Axis(gl[3, 1], xlabel = "", ylabel = L"\beta_x", title = "Beta Functions")
        for (i, data) in enumerate(twi_data)
            lines!(ax_beta_x, floor_length, data[:,1], color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_beta_x)
        
        # 4. Beta Y
        ax_beta_y = Axis(gl[4, 1], xlabel = "", ylabel = L"\beta_y")
        for (i, data) in enumerate(twi_data)
            lines!(ax_beta_y, floor_length, data[:,4], color = beam_colors[i], linewidth = 2)
        end
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
            
            # Use element width
            ele_width = ele.len
            
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
        for (i, data) in enumerate(twi_data)
            lines!(ax_alpha_x, floor_length, data[:,2], color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_alpha_x)
        
        # 7. Alpha Y
        ax_alpha_y = Axis(gl[7, 1], xlabel = "", ylabel = L"\alpha_y")
        for (i, data) in enumerate(twi_data)
            lines!(ax_alpha_y, floor_length, data[:,5], color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_alpha_y)
        
        # 8. Emittance X
        ax_emit_x = Axis(gl[8, 1], xlabel = "", ylabel = L"\epsilon_x", title = "Emittance")
        for (i, data) in enumerate(twi_data)
            lines!(ax_emit_x, floor_length, data[:,3], color = beam_colors[i], linewidth = 2)
        end
        hidexdecorations!(ax_emit_x)
        
        # 9. Emittance Y (bottom plot)
        ax_emit_y = Axis(gl[9, 1], xlabel = "z [m]", ylabel = L"\epsilon_y")
        for (i, data) in enumerate(twi_data)
            lines!(ax_emit_y, floor_length, data[:,6], color = beam_colors[i], linewidth = 2)
        end
        
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