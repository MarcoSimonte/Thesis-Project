#Power spectrum

#Lettura file e conversione dei parametri
using HDF5
using DelimitedFiles
using FFTW
using Statistics
main="/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap="0242"
include(string(main,"read.jl"))
d_conv,d_dm_conv,temp,vx_conv,vy_conv,vz_conv,Npoint=read(snap,main)

#Creazione della griglia
const pc=3.08*10.0^18
const kpc=(10.0^3)*pc
const Nbox=128

#for i in 1:25 # perform the routine at various radii or on box of different size
for i in 1:8
    #Nbox=20+i*10
    d_max,index=findmax(d_conv)
    xc=index[1]
    yc=index[2]
    zc=index[3]

    if (i >= 2)
        xc=xc-(i-1)*10 #moving the center of the different boxes
        yc=yc-(i-1)*10
        zc=zc-(i-1)*10
    end
    rbox=20. *kpc*(i-1.)*10. *sqrt(3.)
    rvir=1.01*1.e3*kpc
    #println(rbox/rvir)
    #Trasformate di Fourier
    main="/home/marco/Scrivania/Tesi/Ammassi_tesi/"
    include(string(main,"FT.jl"))
    vx_k,vy_k,vz_k,vk,dps=fourier_trans(vx_conv,vy_conv,vz_conv,d_conv,xc,yc,zc,Nbox) #fourier transform

    #println(fftfreq(51))
    include(string(main,"ps_profile.jl")) # 3D profile
    vxk_prof,vyk_prof,vzk_prof,dps_prof,vk_prof,k=profile(vx_k,vy_k,vz_k,vk,dps,Nbox)

    include(string(main,"slope.jl")) # compute the slope of the power spectrum
    slope_d,slope_v=slope(k,dps_prof,vxk_prof,vyk_prof,vzk_prof)

    include(string(main,"mean2.jl")) # mean slope
    mean_slope_d,err_slope_d,mean_slope_v,err_slope_v=mean_func(slope_d,slope_v)
    output1=string(main,"Profili_radiali/IT92_2/IT92_2_sloped_radial_profile_8q.txt")
    output2=string(main,"Profili_radiali/IT92_2/IT92_2_slopev_radial_profile_8q.txt")

    out1=open(output1,"a+") do out1
    writedlm( out1, [rbox/rvir mean_slope_d err_slope_d])
    end

    out2=open(output2, "a+") do out2
    writedlm( out2, [rbox/rvir mean_slope_v err_slope_v])
    end
end
#using PyCalls
#using PyPlot

#plot(log10.(k[1:50]),log10.(vzk_prof[1:50]),color="black",linestyle="-",marker=".")
#title("Spettro di potenza della densita R7 (8)",fontsize="18.0")
#xlabel(L"$Log(k)$",fontsize="16.0")
#ylabel(L"$Log(|\rho (k)|^2)$ ",fontsize="16.0")

#plot(k[1:floor(Int64,Nbox/2)],vxk_slope1[1:floor(Int64,Nbox/2)],color="black",marker=".")
#title("Profilo slope della velocità",fontsize="18.0")
#xlabel("nr. bin",fontsize="16.0")
#ylabel("slope ",fontsize="16.0")

#using DelimitedFiles

#output1=string(main,"new_radial_profile_slope_d_5q.txt")
#output2=string(main,"new_radial_profile_slope_v_5q.txt")


#out1=open(output1,"a+") do out1
#writedlm( out1, [r mean_slope_d err_slope_d])
#end

#out2=open(output2, "a+") do out2
#writedlm( out2, [r mean_slope_v err_slope_v])
#end
#exit()
