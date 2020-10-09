#function che fa le traformate di Fourier

function fourier_trans(vx_conv,vy_conv,vz_conv,d_conv,xc,yc,zc,Nbox::Int64)

         N=Nbox/2
         N=floor(Int32, N)
         vx_ft=fft(vx_conv[xc-N:xc+N,yc-N:yc+N,zc-N:zc+N]) # fourier transform of each component of the velocity field
         vy_ft=fft(vy_conv[xc-N:xc+N,yc-N:yc+N,zc-N:zc+N])
         vz_ft=fft(vz_conv[xc-N:xc+N,yc-N:yc+N,zc-N:zc+N])
         d_ft=fft(d_conv[xc-N:xc+N,yc-N:yc+N,zc-N:zc+N])


        # vx_ft_conj=conj(vx_ft)
        # vy_ft_conj=conj(vy_ft)
        # vz_ft_conj=conj(vz_ft)
        # d_ft_conj=conj(d_ft)

         vx_k=zeros(Nbox+1,Nbox+1,Nbox+1)
         vy_k=zeros(Nbox+1,Nbox+1,Nbox+1)
         vz_k=zeros(Nbox+1,Nbox+1,Nbox+1)
         dps=zeros(Nbox+1,Nbox+1,Nbox+1)
         vk=zeros(Nbox+1,Nbox+1,Nbox+1)
         @inbounds for i in 1:Nbox+1
             @inbounds for j in 1:Nbox+1
           @inbounds @simd for k in 1:Nbox+1
                     @fastmath vx_k[i,j,k]=vx_ft[i,j,k]*conj(vx_ft[i,j,k]) #compute vx(k)^2
                     @fastmath vy_k[i,j,k]=vy_ft[i,j,k]*conj(vy_ft[i,j,k])
                     @fastmath vz_k[i,j,k]=vz_ft[i,j,k]*conj(vz_ft[i,j,k])
                     @fastmath vk[i,j,k]=vx_k[i,j,k]+vy_k[i,j,k]+vz_k[i,j,k]
                     @fastmath dps[i,j,k]=d_ft[i,j,k]*conj(d_ft[i,j,k]) # compute density(k)^2
                           end
                       end
                   end

         return(vx_k,vy_k,vz_k,vk,dps)
end
