#calcolo della slope
function slope(k,d,vx,vy,vz)

         slope_d=zeros(0)
         slope_v=zeros(0)
         #slope_vy=zeros(0)
         #slope_vz=zeros(0)
         for i in 1:length(k)-2
               if (log10(k[i]) > 0.55 && log10(k[i]) < 1.25)
                   sloped=(log(d[i+1])-log(d[i]))/(log(k[i+1])-log(k[i])) #slope in each bin of the power spectrum
                   slopevx=(log(vx[i+1])-log(vx[i]))/(log(k[i+1])-log(k[i]))
                   slopevy=(log(vy[i+1])-log(vy[i]))/(log(k[i+1])-log(k[i]))
                   slopevz=(log(vz[i+1])-log(vz[i]))/(log(k[i+1])-log(k[i]))

                   append!(slope_d,sloped)
                   append!(slope_v,[slopevx,slopevy,slopevz])
                   #append!(slope_vy,slopevy)
                   #append!(slope_vz,slopevz)
               end
         end

         return slope_d,slope_v#,slope_vy,slope_vz
end
