#lettura e conversione
function read(snap,root)
         cluster="Ammassi/"
         file1=string(root,cluster,"declust_z",snap)
         file2=string(root,cluster,"declust_v_z",snap)
         file3=string(root,cluster,"deDD",snap,".conv2")

         conv=zeros(5)

         d=h5read(file1,"Density")
         dm=h5read(file1,"Dark_Matter_Density")
         temp=h5read(file1,"Temperature")
         vx=h5read(file2,"x-velocity")
         vy=h5read(file2,"y-velocity")
         vz=h5read(file2,"z-velocity")

         conv=readdlm(file3)

         d_conv=d*conv[4]
         d_dm_conv=dm*conv[4]
         vx_conv=vx*conv[5]
         vy_conv=vy*conv[5]
         vz_conv=vz*conv[5]

         Npoint=length(d[1,1,:])
         #…..Franco: in modo equivalente si può anche usare n=size(d)
         #………………………………………e Npoint=n[3] in questo caso

         return d_conv,d_dm_conv,temp,vx_conv,vy_conv,vz_conv,Npoint, conv[4]
end
