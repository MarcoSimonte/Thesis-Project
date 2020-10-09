#profilo in k della potenza
function profile(d_copy, d_dm, vx, vy, vz, xc, yc, zc,
                 G,kpc, r_max, rmin, Ncell, Ncell_half)#, max_d, min_d)

         #main="/home/marco/Scrivania/Tesi/Ammassi_tesi/"
         radial_bin=25

         kb = 1.28e-16
         mu = 0.62
         mp = 1.67e-24


         radius=zeros(radial_bin)
         for i in 1:radial_bin
             radius[i]=rmin+(i-1.)*r_max/(radial_bin-1.) # uniform grid. Each bin in connected to a radius value
         end

         pk1=zeros(radial_bin)
         pk2=zeros(radial_bin)
         pk3=zeros(radial_bin)
         pk4=zeros(radial_bin)
         N=zeros(radial_bin)
         #pk5=zeros(Nbox+1)
         pk6=zeros(radial_bin)
         #k=zeros(Nbox+1)

         P=zeros(Npoint, Npoint, Npoint)
         v=zeros(Npoint, Npoint, Npoint)
         @inbounds for i in xc-Ncell_half+1:xc+Ncell_half #(xc-N_half:xc+N_half) se voglio una boxcon Ncell centrata in xc,yc,zc)
             @inbounds for j in yc-Ncell_half+1:yc+Ncell_half
           @inbounds @simd for k in zc-Ncell_half+1:zc+Ncell_half
                     @fastmath P[i ,j, k]=d_copy[i, j, k] * kb * temp[i, j, k] / (mu * mp) # compute pressure
                     @fastmath v[i, j, k]=sqrt(vx[i, j, k]^2.0 +vy[i, j, k]^2.0 +vz[i, j, k]^2.0) # compute the modulus of the velocity
                     @fastmath r=Int(floor(sqrt((k -zc -Ncell_half+ 1 -0.5)^2. +(j-yc -Ncell_half +1 -0.5)^2. +(i- xc -Ncell_half+1 -0.5)^2.)))
                               ra=convert(Int64,trunc(r))  # assign radial bin
                               if ra <=1
                                  ra=1
                               end
                               if ra >= radial_bin
                                  ra=radial_bin
                               end

                               pk1[ra]+= d_copy[i,j,k]
                               pk2[ra]+= P[i,j,k]
                               pk3[ra]+= d_dm[i,j,k]
                               pk4[ra]+= 1.

                               #if (d[i, j, k] >= min_d && d[i, j, k] <= max_d)
                                  pk6[ra]+= v[i,j,k]^2.0
                                  N[ra]+=1.
                               #end
                           end
                       end
                   end

         for i in eachindex(pk1)
             pk1[i]=pk1[i]/pk4[i]
             pk2[i]=pk2[i]/pk4[i]
             pk3[i]=pk3[i]/pk4[i]
             #pk4[i]=pk4[i]/pk5[i]
             pk6[i]=sqrt(pk6[i]/N[i])
         end


         mass=zeros(radial_bin)
         g=zeros(radial_bin)
         dr=radius[3]-radius[2]
         for i in eachindex(pk1)
             if(i==1)
                 mass[i]=(4. *pi*(pk1[i]+pk3[i])*radius[i]^3.)/3. # total mass (dark matter + gas)
                 g[i]=G*mass[i]/radius[i]^2.
             else
                 mass[i]=mass[i-1]+(4. *pi*(pk1[i]+pk3[i])*dr*radius[i]^2.) # total mass (daark matter + gas)
                 g[i]=G*mass[i]/radius[i]^2. # "Newton's" gravitational acceleration
             end
         end


         return pk1,pk2,pk3,mass,g,pk6,radius
end
