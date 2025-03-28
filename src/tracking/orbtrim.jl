function pass!(ele::ORBTRIM, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)

    theta_x = ele.theta_x
    theta_y = ele.theta_y
    
    # Calculate kicks if realpara is true
    if ele.realpara
        # C0 is the speed of light in m/s
        # C0 = speed_of_light
        amu = 931494320.0
        # println("charge: ", particles.charge)
        # println("mass (eV/u): ", particles.mass / amu)
        # println("energy (eV/u): ",  (particles.energy / (particles.mass / amu)))
        # println("denom: ",  sqrt(( (particles.energy / (particles.mass / amu)) + amu)^2 - amu^2))
        mass_number = particles.mass/amu
        ecpi = 50.  / (mass_number) * speed_of_light / sqrt(( (particles.energy / mass_number) + amu)^2 - amu^2)
        # println(ecpi)
        theta_x = ele.tm_xkick * ecpi * (particles.charge) / 50.
        theta_y = ele.tm_ykick * ecpi * (particles.charge) / 50.
        
    end
    # Calculate rotation values if needed
    xyrotate_rad = ele.xyrotate * π / 180.0
    
    # Apply kicks to each particle
    for c in 1:num_particles
        if particles.lost_flag[c] == 1
            continue
        end
        
        r6 = @view r_in[(c-1)*6+1:c*6]
        
        # Apply misalignment at entrance
        if !iszero(ele.T1)
            addvv!(r6, ele.T1)
        end
        if !iszero(ele.R1)
            multmv!(r6, ele.R1)
        end
        
        # Apply kicks to momenta
        # println(theta_x[1])
        r6[2] += theta_x
        r6[4] += theta_y
        
        # Apply rotation if needed
        # if ele.xyrotate != 0.0
        #     cos_rot = cos(xyrotate_rad)
        #     sin_rot = sin(xyrotate_rad)
        #     px_temp = r6[2]
        #     py_temp = r6[4]
        #     r6[2] = cos_rot * px_temp - sin_rot * py_temp
        #     r6[4] = sin_rot * px_temp + cos_rot * py_temp
        # end
        
        # Apply misalignment at exit
        if !iszero(ele.R2)
            multmv!(r6, ele.R2)
        end
        if !iszero(ele.T2)
            addvv!(r6, ele.T2)
        end
        
        # Check if particle is lost
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    
    return nothing
end

function pass_TPSA!(ele::ORBTRIM, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # Extract parameters from the element
    theta_x = ele.theta_x
    theta_y = ele.theta_y
    
    # Calculate rotation values if needed
    xyrotate_rad = ele.xyrotate * π / 180.0
    
    # Apply misalignment at entrance
    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r_in, T1)
    end
    if !iszero(R1)
        multmv!(r_in, R1)
    end
    
    # Apply kicks to momenta
    r_in[2] += T(theta_x)
    r_in[4] += T(theta_y)
    
    # Apply rotation if needed
    if ele.xyrotate != 0.0
        cos_rot = cos(xyrotate_rad)
        sin_rot = sin(xyrotate_rad)
        px_temp = r_in[2]
        py_temp = r_in[4]
        r_in[2] = T(cos_rot) * px_temp - T(sin_rot) * py_temp
        r_in[4] = T(sin_rot) * px_temp + T(cos_rot) * py_temp
    end
    
    # Apply misalignment at exit
    if !iszero(R2)
        multmv!(r_in, R2)
    end
    if !iszero(T2)
        addvv!(r_in, T2)
    end   
    
    return r_in
end