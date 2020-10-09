#calcolo di v_turb per il profilo radiale del numero di Ri
function v_turb2(d,vx_conv,vy_conv,vz_conv,kpc,main,Npoint)

    Ncell=15
    N_half=floor(Int32,Ncell/2)
    L=Ncell*20. *kpc
    vx_turb=zeros(Npoint-26,Npoint-26,Npoint-26)
    vy_turb=zeros(Npoint-26,Npoint-26,Npoint-26)
    vz_turb=zeros(Npoint-26,Npoint-26,Npoint-26)
    #v_box= zeros(Npoint-26,Npoint-26,Npoint-26)
@inbounds for i in 26:Npoint-26
@inbounds     for j in 26:Npoint-26
@inbounds @simd   for k in 26:Npoint-26

                      vx_mean, vy_mean, vz_mean=v_mean2(i, j, k, Ncell, N_half, d, vx_conv, vy_conv, vz_conv, Npoint) #evaluating the density/velocity mean on a box wide L
                      vx_turb[i,j,k] = vx_conv[i,j,k]-vx_mean #filtering each velocity component substracting the velocity mean to the original field
                      vy_turb[i,j,k] = vy_conv[i,j,k]-vy_mean
                      vz_turb[i,j,k] = vz_conv[i,j,k]-vz_mean
                      #v = sqrt( vx_conv[i,j,k]^2. +  vy_conv[i,j,k]^2. +  vz_conv[i,j,k]^2.)
                      #v_mean = sqrt( vx_mean^2. + vy_mean^2. + vz_mean^2.)
                      #v_box[i,j,k]= log(v / v_mean) #evaluating the logarithmic velocity perturbations
                  end
              end
          end
    return(vx_turb, vy_turb, vz_turb, L)
end
