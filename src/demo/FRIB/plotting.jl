# Run in a new REPL for each use

include("/Users/Keeney/JuTrack.jl/src/JuTrack.jl")
using .JuTrack

using GLMakie
using Serialization
using Observables
using Gtk

# Lots of this is derived from the Makie documentation at https://docs.makie.org/stable/

# GUI Function
function lattice_gui()
    #############################################
    #                ~ Lattice ~                #
    #############################################

    # Function to plot the labels of each different type within the lattice
    function labels(lengths_list, label, axis, lens)
        texts = []
        for i in eachindex(lengths_list)
            text = text!(axis, label[i], position=((lengths_list[i] - lens[i]/2), 0), space = :data, fontsize = 10)
            push!(texts, text)
        end
        return texts
    end

    # Function to allow the user to turn the labels on or off
    function toggle_labels(labels_list, toggles)
        for label in labels_list
            obs_visible = Observable(label.visible)
            connect!(obs_visible[], toggles[1].active)
        end
    end

    # Gets the different values from the twissring function: beta functions, dispersion, and tune values
    function ring_functions(ring)
        twi = twissring(ring, 0.0, 1)

        twi_betax = []
        for i in eachindex(twi)
            append!(twi_betax, twi[i].betax)
        end

        twi_betay = []
        for i in eachindex(twi)
            append!(twi_betay, twi[i].betay)
        end

        twi_disp = []
        for i in eachindex(twi)
            append!(twi_disp, twi[i].dx * 60)
        end

        tunex = twi[end].dmux/2/pi
        tuney = twi[end].dmuy/2/pi

        return twi_betax, twi_betay, twi_disp, tunex, tuney
    end

    function plot_ring_functions(ring, f, ax, s)
        twi_betax, twi_betay, twi_disp, tunex, tuney = ring_functions(ring)

        toggles = [Toggle(f, active = active, buttoncolor = RGBAf(0,0,0,0.8), framecolor_active = RGBAf(0,0,0,0.2)) for active in [true, false, false]]
        labels = [Label(f, lift(x -> x ? "$l function" : "$l function", t.active))
        for (t, l) in zip(toggles, ["Beta x", "Beta y", "60*Dispersion"])]
        f[3, 1] = grid!(hcat(toggles, labels), tellheight = false)

        # convert function help from https://stackoverflow.com/questions/35482527/how-do-i-change-the-data-type-of-a-julia-array-from-any-to-float64
        bx = lines!(ax, s, convert(Array{Float32,1}, twi_betax))
        by = lines!(ax, s, convert(Array{Float32,1}, twi_betay))
        d = lines!(ax, s, convert(Array{Float32,1}, twi_disp), color = :firebrick)

        Label(f[4, 1], "Tune x = $tunex")
        Label(f[5, 1], "Tune y = $tuney")

        connect!(bx.visible, toggles[1].active)
        connect!(by.visible, toggles[2].active)
        connect!(d.visible, toggles[3].active)
    end

    f = Figure(size = (850, 650))
    f1 = Figure(size = (850, 650))
    f2 = Figure(size = (850, 650))
    ax = Axis(f1[1:8, 2], limits = (nothing, nothing, -1.0, 10.0), xlabel = "s (m)", ylabel = "Beta", title = "Lattice")

    # Creating the lattice
    function lattice(file, f, ax)
        # Load file
        Lat = deserialize(file)

        # Turns the matrix into one long list (in order)
        flattened_ring = vcat(Lat...)

        #Initialize the graph
        s = spos(flattened_ring)

        # Calculates beta functions + dispersion if applicable
        functions_button = Button(f[1, 1], label = "Calculate \n Beta and Dispersion \n Functions", buttoncolor_hover = RGBAf(0,0,0,0.15))
        on(functions_button.clicks) do _
            plot_ring_functions(Lat, f, ax, s)
        end

        # Different colors for the lattice
        colors = [:white, :firebrick3, :skyblue, :green, :yellow, :orange, :purple, :gray50, :mediumvioletred, :chocolate4, :teal]

        # Get each separate length value in order
        lens = []
        for i in eachindex(s)
            if i == 1
                append!(lens, s[i])
                continue
            else
                val = s[i] - s[i-1]
                append!(lens, val)
                continue
            end
        end

        # Set up the toggles to turn the labels on/off
        toggles = [Toggle(f, active = active, buttoncolor = RGBAf(0,0,0,0.8), framecolor_active = RGBAf(0,0,0,0.2)) for active in [false]]
        names = [Label(f, lift(x -> x ? "$l on" : "$l off", t.active))
        for (t, l) in zip(toggles, ["Labels"])]
        f[6, 1] = grid!(hcat(toggles, names), tellheight = false)

        # Find the positons of each different elements within the "Line" list
        drifts = findelem(flattened_ring, :eletype, "DRIFT")
        quadrupoles = findelem(flattened_ring, :eletype, "KQUAD")
        sextupoles = findelem(flattened_ring, :eletype, "KSEXT")
        octupoles = findelem(flattened_ring, :eletype, "KOCT")
        multipoles = findelem(flattened_ring, :eletype, "thinMULTIPOLE")
        solenoids = findelem(flattened_ring, :eletype, "SOLENOID")
        rfcas = findelem(flattened_ring, :eletype, "RFCA")
        sbends = findelem(flattened_ring, :eletype, "SBEND")
        rbends = findelem(flattened_ring, :eletype, "RBEND")
        hkickers = findelem(flattened_ring, :eletype, "HKICKER")
        vkickers = findelem(flattened_ring, :eletype, "VKICKER")
        markers = findelem(flattened_ring, :eletype, "MARKER")

        # Separate quadrupoles based on their k1 value
        quads_pos = []
        quads_neg = []
        other = []
        global count = 1
        for val in flattened_ring
            if val.eletype == "KQUAD"
                if val.k1 >= 0
                    push!(quads_pos, count)
                else
                    push!(quads_neg, count)
                end
            else
                push!(other, count)
            end
            global count += 1
        end

        # Create the plot for the lattice
        for i in eachindex(s)
            if i == 1 # Plots the first element of the lattice
                # Equations to plot the shapes
                rect1_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 1), ((3s[i] - lens[i])/2, 1), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                rect2_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, -1), ((3s[i] - lens[i])/2, -1), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                small_rect_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                triangle1_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 1), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                triangle2_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((3s[i] - lens[i])/2, 1), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                triangle3_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i])/2, 0.9), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                pentagon1_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.8), ((s[i])/2, 0.6), ((3s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                pentagon2_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.8), ((s[i])/2, 1), ((3s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                pentagon3_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.9), ((3s[i] - lens[i])/2, 0.9), ((s[i])/2, 0.45), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                pentagon4_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i])/2, 0.45), ((s[i] - lens[i])/2, 0.9), ((3s[i] - lens[i])/2, 0.9), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                trapezoid1_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.65), ((3s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                trapezoid2_first_coords = Point2f[((s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0.8), ((3s[i] - lens[i])/2, 0.65), ((3s[i] - lens[i])/2, 0), ((s[i] - lens[i])/2, 0)]
                
                # Plotting based on which list the index is in
                if i in drifts
                    poly!(ax, rect1_first_coords, color = colors[1], alpha = 0)
                elseif i in quads_pos
                    poly!(ax, rect1_first_coords, color = colors[2], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in quads_neg
                    poly!(ax, rect2_first_coords, color = colors[2], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in sextupoles
                    poly!(ax, trapezoid1_first_coords, color = colors[3], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in octupoles
                    poly!(ax, trapezoid2_first_coords, color = colors[4], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in multipoles
                    poly!(ax, small_rect_first_coords, color = colors[5], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in solenoids
                    poly!(ax, triangle2_first_coords, color = colors[6], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in rfcas
                    poly!(ax, triangle3_first_coords, color = colors[7], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in sbends
                    poly!(ax, pentagon3_first_coords, color = colors[8], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in rbends
                    poly!(ax, pentagon4_first_coords, color = colors[8], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in hkickers
                    poly!(ax, pentagon1_first_coords, color = colors[9], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in vkickers
                    poly!(ax, pentagon2_first_coords, color = colors[10], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in markers
                    poly!(ax, triangle1_first_coords, color = colors[11], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                end
            else # Plots the rest of the elements in the lattice
                # More shapes
                rect1_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 1), ((s[i] - lens[i]/2) + (lens[i]/2), 1), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                rect2_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), -1), ((s[i] - lens[i]/2) + (lens[i]/2), -1), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                small_rect_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                triangle1_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 1), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                triangle2_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) + (lens[i]/2), 1), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                triangle3_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2), 0.9), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                pentagon1_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.8), ((s[i] - lens[i]/2), 0.6), ((s[i] - lens[i]/2) + (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                pentagon2_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.8), ((s[i] - lens[i]/2), 1), ((s[i] - lens[i]/2) + (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                pentagon3_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.9), ((s[i] - lens[i]/2) + (lens[i]/2), 0.9), ((s[i] - lens[i]/2), 0.45), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                pentagon4_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2), 0.45), ((s[i] - lens[i]/2) - (lens[i]/2), 0.9), ((s[i] - lens[i]/2) + (lens[i]/2), 0.9), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                trapezoid1_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.65), ((s[i] - lens[i]/2) + (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]
                trapezoid2_other_coords = Point2f[((s[i] - lens[i]/2) - (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0.8), ((s[i] - lens[i]/2) + (lens[i]/2), 0.65), ((s[i] - lens[i]/2) + (lens[i]/2), 0), ((s[i] - lens[i]/2) - (lens[i]/2), 0)]

                if i in drifts # Help from ChatGPT setting the colors based on the type of element it is
                    poly!(ax, rect1_other_coords, color = colors[1], alpha = 0)
                elseif i in quads_pos
                    poly!(ax, rect1_other_coords, color = colors[2], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in quads_neg
                    poly!(ax, rect2_other_coords, color = colors[2], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in sextupoles
                    poly!(ax, trapezoid1_other_coords, color = colors[3], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in octupoles
                    poly!(ax, trapezoid2_other_coords, color = colors[4], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in multipoles
                    poly!(ax, small_rect_other_coords, color = colors[5], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in solenoids
                    poly!(ax, triangle2_other_coords, color = colors[6], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in rfcas
                    poly!(ax, triangle3_other_coords, color = colors[7], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in sbends
                    poly!(ax, pentagon3_other_coords, color = colors[8], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in rbends
                    poly!(ax, pentagon4_other_coords, color = colors[8], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in hkickers
                    poly!(ax, pentagon1_other_coords, color = colors[9], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in vkickers
                    poly!(ax, pentagon2_other_coords, color = colors[10], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                elseif i in markers
                    poly!(ax, triangle1_other_coords, color = colors[11], strokecolor = :black, strokewidth = 0.7, alpha = 0.6)
                end
            end

            # Plot each center value, and allow the points to be turned on/off by the toggle
            scatter = scatter!(ax, (s[i] - (lens[i]/2)), 0, color = :black) 
            connect!(scatter.visible, toggles[1].active) # Make scatter invisible if toggle is off
        
        end

        # Putting each user-defined name into a list in order
        names = []
        for val in flattened_ring
            push!(names, val.name)
        end
        names = convert(Vector{String}, names)

        # Plotting the labels in the correct positions
        labels_names = labels(s, names, ax, lens)

        # Toggle can turn on/off each different label
        toggle_labels(labels_names, toggles)

        # Creating the Legend
        elem_quad = [PolyElement(polycolor = (colors[2], 0.6), strokecolor = :black, strokewidth = 0.7)]
        elem_sext = [PolyElement(polycolor = (colors[3], 0.6), strokecolor = :black, strokewidth = 0.7,
                    points = Point2f[(0, 0), (0, 0.65), (1, 0.8), (1, 0), (0, 0)])]
        elem_oct = [PolyElement(polycolor = (colors[4], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 0.8), (1, 0.65), (1, 0), (0, 0)])]
        elem_mult = [PolyElement(polycolor = (colors[5], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 0.8), (1, 0.8), (1, 0), (0, 0)])]
        elem_sol = [PolyElement(polycolor = (colors[6], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (1, 1), (1, 0), (0, 0)])]
        elem_rfca = [PolyElement(polycolor = (colors[7], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0.5, 0.9), (1, 0), (0, 0)])]
        elem_bend = [PolyElement(polycolor = (colors[8], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 0.9), (1, 0.9), (0.5, 0.45), (1, 0), (0, 0)])]
        elem_hkick = [PolyElement(polycolor = (colors[9], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 0.8), (0.5, 0.6), (1, 0.8), (1, 0), (0, 0)])]
        elem_vkick = [PolyElement(polycolor = (colors[10], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 0.8), (0.5, 1), (1, 0.8), (1, 0), (0, 0)])]
        elem_marker = [PolyElement(polycolor = (colors[11], 0.6), strokecolor = :black, strokewidth = 0.7, 
                    points = Point2f[(0, 0), (0, 1), (1, 0), (0, 0)])]

        Legend(f1[1:8, 3], 
        [elem_quad, elem_sext, elem_oct, elem_mult, elem_sol, elem_rfca, elem_bend, elem_hkick, elem_vkick, elem_marker], 
        ["Quadrupole", "Sextupole", "Octupole", "Thin Multipole", "Solenoid", "RFCA", "Bend", "HKicker", "VKicker", "Marker"], 
        "Legend")

        # Displaying the lattice
        display(f)
    end

    #############################################
    #             ~ Beam Tracking ~             #
    #############################################
    # Function so that the value placed in the textbox will change in the calculation
    function change_tb(tb, observable, type)
        on(tb.stored_string) do s
            observable[] = parse(type, s)
        end
    end

    # Function so that the textbox settings don't clutter the code as much
    function shorter_textbox(layout, placeholder, validator)
        textbox = Textbox(layout, placeholder = placeholder, validator = validator, tellwidth = false, 
        bordercolor = :black, bordercolor_focused = :black, cornerradius = 2, width = 70, boxcolor_focused_invalid = RGBAf(1,0,0,0.15), 
        bordercolor_focused_invalid = RGBAf(0,0,0,1))
        return textbox
    end

    # Values in the function declaration are the default values in the textboxes
    function beam_plot(lattice, f, nmacro = 1, nreal = 1, e0 = 10e9, ϕs = 10.0, v = 2.9e7, betax = 0.50, betay = 0.50, harmonic_number = 7560.0, freq = 591e6, αc = 3.5, emitx = 2e-9, emity = 2e-9, emitz = 2e-4, nturns = 5)
        Lat = deserialize(lattice)
        flattened_ring = vcat(Lat...)

        # Defining layouts for positioning
        layout_left = GridLayout()

        ax0 = Axis(f[1:18, 1:2], tellwidth = false, height = 100)
        layout_left[1:18, 1:2] = ax0
        f.layout[1:18, 1:2] = layout_left
        hidedecorations!(ax0)
        hidespines!(ax0)

        # Adding the textboxes and the labels on the outside
        label_header = Label(layout_left[2, 1:2], "Variables", fontsize = 20)

        tb_nmacro = shorter_textbox(layout_left[3, 2], "1", Int)
        label_nmacro = Label(layout_left[3, 1], "N Macro =", justification = :right, width = 110, halign = :left)

        tb_nreal = shorter_textbox(layout_left[4, 2], "1", Float64)
        label_nreal = Label(layout_left[4, 1], "N Real =", justification = :right, width = 110)

        tb_e0 = shorter_textbox(layout_left[5, 2], "10e9", Float64)
        label_e0 = Label(layout_left[5, 1], "E0 (eV) =", justification = :right, width = 110)

        tb_ϕs = shorter_textbox(layout_left[6, 2], "10.0", Float64)
        label_ϕs = Label(layout_left[6, 1], "ϕs (deg) =", justification = :right, width = 110)

        tb_v = shorter_textbox(layout_left[7, 2], "2.9e7", Float64)
        label_v = Label(layout_left[7, 1], "Voltage (Volt) =", justification = :right, width = 110)

        tb_betax = shorter_textbox(layout_left[8, 2], "0.50", Float64)
        label_betax = Label(layout_left[8, 1], "Beta x (m) =", justification = :right, width = 110)

        tb_betay = shorter_textbox(layout_left[9, 2], "0.50", Float64)
        label_betay = Label(layout_left[9, 1], "Beta y (m) =", justification = :right, width = 110)

        tb_hn = shorter_textbox(layout_left[10, 2], "7560.0", Float64)
        label_hn = Label(layout_left[10, 1], "Harmonic Number =", justification = :right, width = 110)

        tb_freq = shorter_textbox(layout_left[11, 2], "591e6", Float64)
        label_freq = Label(layout_left[11, 1], "Frequency (Hz) =", justification = :right, width = 110)

        tb_αc = shorter_textbox(layout_left[12, 2], "3.5", Float64)
        label_αc = Label(layout_left[12, 1], "αc =", justification = :right, width = 110)

        tb_emitx = shorter_textbox(layout_left[13, 2], "2e-9", Float64)
        label_emitx = Label(layout_left[13, 1], "Emit x (m*rad) =", justification = :right, width = 110)

        tb_emity = shorter_textbox(layout_left[14, 2], "2e-9", Float64)
        label_emity = Label(layout_left[14, 1], "Emit y (m*rad) =", justification = :right, width = 110)

        tb_emitz = shorter_textbox(layout_left[15, 2], "2e-4", Float64)
        label_emitz = Label(layout_left[15, 1], "Emit z (s*eV) =", justification = :right, width = 110)

        tb_nturns = shorter_textbox(layout_left[16, 2], "5", Int)
        label_nturns = Label(layout_left[16, 1], "N Turns =", justification = :right, width = 110)

        before_label = Label(f[1, 3:8], "Electron Distribution Before", fontsize = 20)
        after_label = Label(f[8, 3:8], "Electron Distribution After", fontsize = 20)

        # Define observables (initial values for the variables)
        nmacro_obs = Observable(nmacro)
        nreal_obs = Observable(nreal)
        e0_obs = Observable(e0)
        ϕs_obs = Observable(ϕs)
        volt_obs = Observable(v)
        betax_obs = Observable(betax)
        betay_obs = Observable(betay)
        harmonic_obs = Observable(harmonic_number)
        freq_obs = Observable(freq)
        αc_obs = Observable(αc)
        emitx_obs = Observable(emitx)
        emity_obs = Observable(emity)
        emitz_obs = Observable(emitz)
        nturns_obs = Observable(nturns)

        ebeam = Observable(zeros(Float64, nmacro_obs[], 6))
        
        # Allows the input in the textbox to change the value of the variable
        change_tb(tb_nmacro, nmacro_obs, Int)
        change_tb(tb_nreal, nreal_obs, Float64)
        change_tb(tb_e0, e0_obs, Float64)
        change_tb(tb_ϕs, ϕs_obs, Float64)
        change_tb(tb_v, volt_obs, Float64)
        change_tb(tb_betax, betax_obs, Float64)
        change_tb(tb_betay, betay_obs, Float64)
        change_tb(tb_hn, harmonic_obs, Float64)
        change_tb(tb_freq, freq_obs, Float64)
        change_tb(tb_αc, αc_obs, Float64)
        change_tb(tb_emitx, emitx_obs, Float64)
        change_tb(tb_emity, emity_obs, Float64)
        change_tb(tb_emitz, emitz_obs, Float64)
        change_tb(tb_nturns, nturns_obs, Int)

        # Initialize axes
        ax1 = Axis(f[2:7, 3:4], tellwidth = false)
        ax2 = Axis(f[2:7, 5:6], tellwidth = false)
        ax3 = Axis(f[2:7, 7:8], tellwidth = false)

        ax4 = Axis(f[9:15, 3:4], tellwidth = false)
        ax5 = Axis(f[9:15, 5:6], tellwidth = false)
        ax6 = Axis(f[9:15, 7:8], tellwidth = false)

        # Setting initial axis labels
        ax1.title = "x"
        ax2.title = "y"
        ax3.title = "z"
        ax4.title = "x"
        ax5.title = "y"
        ax6.title = "z"

        ax1.xlabel = "x position"
        ax1.ylabel = "px"
        ax2.xlabel = "y position"
        ax2.ylabel = "py"
        ax3.xlabel = "z position"
        ax3.ylabel = "pz"
        ax4.xlabel = "x position"
        ax4.ylabel = "px"
        ax5.xlabel = "y position"
        ax5.ylabel = "py"
        ax6.xlabel = "z position"
        ax6.ylabel = "pz"

        ax1.ylabelrotation = 0
        ax2.ylabelrotation = 0
        ax3.ylabelrotation = 0
        ax4.ylabelrotation = 0
        ax5.ylabelrotation = 0
        ax6.ylabelrotation = 0

        # Makes the graphs easier to change
        function update_plot_before(beam)

            scatter!(ax1, beam[:, 1], beam[:, 2], color = :firebrick2)
            scatter!(ax2, beam[:, 3], beam[:, 4], color = :firebrick2)
            scatter!(ax3, beam[:, 5], beam[:, 6], color = :firebrick2)

            display(f)
        end

        function update_plot_after(beam)

            scatter!(ax4, beam[:, 1], beam[:, 2], color = :firebrick2)
            scatter!(ax5, beam[:, 3], beam[:, 4], color = :firebrick2)
            scatter!(ax6, beam[:, 5], beam[:, 6], color = :firebrick2)

            display(f)
        end

        button = Button(layout_left[17, 2], label = "Calculate", buttoncolor_hover = RGBAf(0,0,0,0.15), cornerradius = 2)
        parallel_button = Button(layout_left[17, 1], label = "Calculate Using \n Parallel \n Computing", buttoncolor_hover = RGBAf(0,0,0,0.15), cornerradius = 2)

        # Update the calculation for when you hit the button
        on(button.clicks) do _
            # Clear the axes whenever the "Calculate" button is clicked to allow for new results
            for ax in (ax1, ax2, ax3, ax4, ax5, ax6)
                empty!(ax)
            end

            # Update the graphs based on the given values
            mainRFe = AccelCavity(freq_obs[], volt_obs[], harmonic_obs[], π-ϕs_obs[]*π/180.0)
            lmap = LongitudinalRFMap(αc_obs[], mainRFe)
            opIPe = optics4DUC(betax_obs[], 0.0, betay_obs[], 0.0)
            beam = Beam(zeros(nmacro_obs[], 6), np = nreal_obs[], energy = e0_obs[], emittance=[emitx_obs[], emity_obs[], emitz_obs[]])
            initilize_6DGaussiandist!(beam, opIPe, lmap)
            get_emittance!(beam)
            emit_before = beam.emittance
            Box(f[16, 3:8], strokecolor = :white)
            Label(f[16, 3:8], "Emittance before: $emit_before")

            beam_data_before = fetch(beam.r)
            ebeam[] = beam_data_before
            update_plot_before(fetch(beam_data_before))

            ringpass!(flattened_ring, beam, nturns_obs[])
            beam_data_after = fetch(beam.r)
            ebeam[] = beam_data_after
            get_emittance!(beam)
            emit_after = beam.emittance
            Box(f[17, 3:8], strokecolor = :white)
            Label(f[17, 3:8], "Emittance after: $emit_after")
            update_plot_after(fetch(beam_data_after))
        end

        # Same thing as above, but with parallel computing
        on(parallel_button.clicks) do _
            # Clear the axes whenever the "Calculate" button is clicked to allow for new results
            for ax in (ax1, ax2, ax3, ax4, ax5, ax6)
                empty!(ax)
            end

            # Update the graphs based on the given values
            mainRFe = AccelCavity(freq_obs[], volt_obs[], harmonic_obs[], π-ϕs_obs[]*π/180.0)
            lmap = LongitudinalRFMap(αc_obs[], mainRFe)
            opIPe = optics4DUC(betax_obs[], 0.0, betay_obs[], 0.0)
            beam = Beam(zeros(nmacro_obs[], 6), np = nreal_obs[], energy = e0_obs[], emittance=[emitx_obs[], emity_obs[], emitz_obs[]])
            initilize_6DGaussiandist!(beam, opIPe, lmap)
            get_emittance!(beam)
            emit_before = beam.emittance
            Box(f[16, 3:8], strokecolor = :white)
            Label(f[16, 3:8], "Emittance before: $emit_before")

            beam_data_before = fetch(beam.r)
            ebeam[] = beam_data_before
            update_plot_before(fetch(beam_data_before))

            pringpass!(flattened_ring, beam, nturns_obs[])
            beam_data_after = fetch(beam.r)
            ebeam[] = beam_data_after
            get_emittance!(beam)
            emit_after = beam.emittance
            Box(f[17, 3:8], strokecolor = :white)
            Label(f[17, 3:8], "Emittance after: $emit_after")
            update_plot_after(fetch(beam_data_after))
        end

        update_plot_before(zeros(nmacro_obs[], 6))
    end

    # Information + Layout for the first screen
    Box(f[1:11, 1:5], color = :white, tellwidth = false)
    Box(f[1:11, 6:10], color = :white, tellwidth = false)
    Label(f[2, 3], "Lattice Display", fontsize = 35)
    Label(f[2, 8], "Beam Tracking", fontsize = 35)
    Label(f[5, 1:5], "Displays the layout of the lattice", fontsize = 18)
    Label(f[6, 1:5], "Graphs Beta and Dispersion functions \n when applicable", fontsize = 18)
    Label(f[7, 1:5], "Calculates the tune values", fontsize = 18)
    Label(f[4, 6:10], "Calculates particle distribution through \n a ring, and displays the distribution", fontsize = 18)
    Label(f[5, 6:10], "Ability to change variables", fontsize = 18)
    Label(f[6, 6:10], "Displays Emittance values", fontsize = 18)
    Label(f[7, 6:10], "Allows Parallel Computing", fontsize = 18)
    Label(f[8, 6:10], "(Note: Click 'Enter' after inputting a value \n into a textbox to register it properly)", fontsize = 15)


    # Buttons to open a file explorer to allow the user to input a .jls lattice file
    button_import_lat = Button(f[10, 3], label = "Import Lattice for \n Display", buttoncolor_hover = RGBAf(0,0,0,0.15))
    button_import_track = Button(f[10, 8], label = "Import Lattice for \n Beam Tracking", buttoncolor_hover = RGBAf(0,0,0,0.15))

    warning_label = Label(f[12, 5:6], "Larger lattices may take a long time to load,\n so please be patient.")

    on(button_import_lat.clicks) do _
        empty!(ax)
        file = open_dialog_native("Choose the input Julia file", GtkNullContainer(), ("*.jls",)) # function from https://stackoverflow.com/questions/73959829/taking-an-excel-file-as-a-input-data-from-user-in-julia#:~:text=You%20can%20use%20open_dialog_native%20from%20Gtk.jl.%20julia%3E%20open_dialog_native%28%22Choose,the%20chosen%20file%27s%20full%20path%20as%20a%20string.
        "/path/to/myfile.jls"
        lattice(file, f1, ax)
    end

    on(button_import_track.clicks) do _
        empty!(ax)
        file = open_dialog_native("Choose the input Julia file", GtkNullContainer(), ("*.jls",))
        "/path/to/myfile.jls"
        beam_plot(file, f2)
    end

    display(f)
end

lattice_gui()