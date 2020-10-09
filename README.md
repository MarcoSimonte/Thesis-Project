# Thesis-Project

In this repository I will store my Julia's script which I used to analyze the data set obtained by cosmological simulations developed with ENZO code.

1. Richardson number and density fluctuations

I tried to evaluate the Richardson number and the strength of the density fluctuations in a simulated ICM. I evaluated these quantity in various boxes at different 
radii. 

Main: "density_fluct_2.jl"
Functions: "read.jl" reads and HDF5 file. 
           "v_turb.jl" , "v_mean.jl" filter the density and velocity field in order to study the properties of the turbulent motions. To achieve this I substracted
           the mean, evaluated in a 300 kpc wide, to the original field. 
           "pdf.jl" computes a probability distribution function.
           "fluct_out.jl" rules out the 5% denser density fluctuations.
           "radial_profile.jl" performs a 3D radiale profile.
           "BV_freq.jl" estimates the Brunt-Vaisala frequency.
           "Ri_box.jl" estimates the Richardson number.
           
2. Power spectrum radiale profile

In this routine I generated velocity and density power spectra. Besides, I estimated the mean slope of this spectra in order to study their radial profile. 

Main: "Power_spectrum.jl" 
Functions: "read.jl" reads and HDF5 file.
           "FT.jl" performs Fourier transform.
           "ps_profile.jl" computes a 3D radial profile in order to create a power spectrum.
           "slope.jl" evaluates the slope of the power spectrum in each bin.
           "mean2.jl" evaluates the mean slope of the power spectrum.
           
