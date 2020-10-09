function d_out(n, bin, main, Ntot, d_max)

        Npoint=length(bin)
        d_max, index=findmax(bin)



        dist=zeros(Npoint)
        for i in 1:Npoint
            dist[i]=abs(bin[i]-d_max)
        end





@inbounds for i in 1:Npoint   #sorting the bin of the pdf considering the distance of the bin from the mode
              imin=i
              distmin=dist[i]
              Nmin=n[i]
              binmin=bin[i]
@inbounds @simd  for j in i:Npoint
                     if(distmin >= dist[j])
                        Nmin=n[j]
                        distmin=dist[j]
                        binmin=bin[j]
                        imin=j
                     end
                 end
                dist[imin]=dist[i]
                n[imin]=n[i]
                bin[imin]=bin[i]
                dist[i]=distmin
                n[i]=Nmin
                bin[i]=binmin
        end


        output2=string(main,"fluttuazioni2.txt")
        out2=open(output2,"w") do out2
        writedlm( out2, [dist n])
        end

        cumulative=0.
        i99=1
        i97=1
        i95=1
        for i in 1:Npoint  # find the 10, 5, 3 % denser fluctuations
                #println(n[i], " ", bin[i])
                cumulative+= n[i]/Ntot
                if (cumulative >= 10. /100.)
                        i99=i
                        break
                end
        end

        cumulative=0.
        for i in 1:Npoint
                cumulative+= n[i]/Ntot
                if (cumulative >= 3. /100.)
                        i97=i
                        break
                end
        end

        cumulative=0.
        for i in 1:Npoint
                cumulative+= n[i]/Ntot
                if (cumulative >= 5. /100.)
                        i95=i
                        break
                end
        end


        return(n,bin,i99,i97, i95)


end
