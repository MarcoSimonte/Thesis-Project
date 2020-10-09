#profilo in k della potenza
function profile(vx_k,vy_k,vz_k,vk,dps,Nbox)

         pk1=zeros(Nbox+1)
         pk2=zeros(Nbox+1)
         pk3=zeros(Nbox+1)
         pk4=zeros(Nbox+1)
         pk4=zeros(Nbox+1)
         pk5=zeros(Nbox+1)
         pk6=zeros(Nbox+1)
         #k=zeros(Nbox+1)

         @inbounds for i in 1:Nbox+1
             @inbounds for j in 1:Nbox+1
           @inbounds @simd for k in 1:Nbox+1
                     @fastmath r=Int(floor(sqrt((k -1. -0.5)^2. +(j-1. -0.5)^2. +(i-1. -0.5)^2.)))
                               ra=convert(Int64,trunc(r))  # assign radial bin
                               if ra <=1
                                  ra=1
                               end
                               if ra >= Nbox+1
                                  ra=Nbox+1
                               end
                               pk1[ra]+= vx_k[i,j,k]
                               pk2[ra]+= vy_k[i,j,k]
                               pk3[ra]+= vz_k[i,j,k]
                               pk4[ra]+= dps[i,j,k]
                               pk5[ra]+= vk[i,j,k]
                               pk6[ra]+= 1
                           end
                       end
                   end
         k=zeros(Nbox+1)
         for i in eachindex(pk1)
             pk1[i]=pk1[i]/pk6[i]
             pk2[i]=pk2[i]/pk6[i]
             pk3[i]=pk3[i]/pk6[i]
             pk4[i]=pk4[i]/pk6[i]
             pk5[i]=pk5[i]/pk6[i]
             k[i]=1.0*i
         end

         return pk1,pk2,pk3,pk4,pk5,k
end
