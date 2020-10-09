using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics

main = "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap = "0199"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

rvir=0.78*1.e3*kpc

include(string(main, "read.jl"))
d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap, main)
Nhalf = floor(Int64, Npoint / 2)

vx_copy = vx
vy_copy = vy
vz_copy = vz
d_max, index = findmax(d)
xc = index[1]
yc = index[2]
zc = index[3]

include(string(main, "v_turb2.jl")) #Filtering the original (density and velocity) field with
include(string(main, "v_mean2.jl")) #fixed scale L
vx, vy, vz, L = v_turb2(d, vx, vy, vz, kpc, main, Npoint)
include(string(main, "v_turb.jl"))
include(string(main, "v_mean.jl"))
vx, vy, vz, L, d_turb, out = v_turb(d, vx, vy, vz, kpc, main, Npoint)


Ncell=32
Ncell_half=floor(Int32, Ncell/2.)

L=300. *kpc
out=26
d=d*conv_d

for i in 1:50 #perform the routine on various boxes at different radii

xc_2= xc -(i-1)*4 #finding the  coordinate of thecenter of the box
yc_2= yc -(i-1)*4
zc_2= zc -(i-1)*4

rad= 20. *kpc* sqrt( (xc_2-xc)^2.0 + (yc_2-yc)^2.0 + (zc_2 - zc)^2.0) #distance of the center of the considered box from the cluster center

if(rad/rvir > 1.)
        exit()
end


include(string(main, "pdf.jl")) #computing the probability distribution function
n, bin, d_max = pdf(d_turb[xc_2-Ncell_half+1:xc_2+Ncell_half,
                           yc_2-Ncell_half+1:yc_2+Ncell_half,
                           zc_2-Ncell_half+1:zc_2+Ncell_half], 500)


Ntot=sum(n)
include(string(main, "fluct_out.jl")) #ruling out the denser fluctuations
n, bin, i99, i97, i95=d_out(n,bin, main, Ntot, d_max)

max_d,ind1=findmax((bin[i95:length(bin)]))
min_d,ind2=findmin((bin[i95:length(bin)]))

dev= var((bin[i95+1:length(bin)])) # variance of the pdf
#include(string(main, "pdf.jl"))

#var_v=var(v_box[xc_2-Ncell_half+1:xc_2+Ncell_half,
#                           yc_2-Ncell_half+1:yc_2+Ncell_half,
#                           zc_2-Ncell_half+1:zc_2+Ncell_half])

include(string(main, "radial_profile.jl"))
include(string(main, "BV_freq.jl"))
include(string(main, "Ri_box.jl"))
rmin= 20. *kpc* sqrt((xc_2-Ncell_half+1 - xc)^2. + (yc_2-Ncell_half+1 - yc)^2. + (zc_2-Ncell_half+1 - zc)^2.)
rmax= 20. *kpc* sqrt((xc_2+Ncell_half - xc)^2. + (yc_2+Ncell_half - yc)^2. + (zc_2+Ncell_half - zc)^2.)


Ri, err_Ri=Ri_box(d, d_turb, d_dm, vx, vy, vz,                   # estimating the Richardson number
                  xc_2, yc_2, zc_2,G, kpc, rmax, rmin, Ncell,
                  Ncell_half, max_d, min_d, L)

output1=string(main,"Rich_number_profile/IT90_4/IT90_4_dfluct_vs_Ri_box_1q.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [rad/rvir  dev Ri err_Ri])
end

end
