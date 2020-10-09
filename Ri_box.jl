function Ri_box(d_conv, d_turb, d_dm_conv, vx, vy, vz,
                xc, yc, zc, G, kpc, rmax, rmin, Ncell,
                Ncell_half, max_d, min_d, L)


    d_prof, P_prof, d_dm_prof, mass,g, v_prof, r =profile(d_conv, d_turb, d_dm_conv, vx, vy, vz,  #compute a 3D radial profile
                                                           xc, yc, zc, G, kpc, rmax, rmin, Ncell,
                                                           Ncell_half, max_d, min_d)


    for i = 1:length(v_prof)
        v_prof[i] = abs(v_prof[i])
    end


    N_square = BV(P_prof, d_prof, r, g, gamma)

    Rich_number = zeros(length(N_square)+1)

    @inbounds for i = 6:length(N_square)
        @fastmath Rich_number[i] = abs(N_square[i]) / ((v_prof[i] / L)^2.0) #definition of the richardson number
                  println(i, Rich_number[i])
        #@fastmath Rich_number2[i] = abs(N_square[i]) / (1.e8 / (1.e3 * kpc))^2.0
    end



    return mean(Rich_number), std(Rich_number)

end
