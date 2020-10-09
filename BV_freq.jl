#calcolo frequenza BV

function BV(P,d,r,g,gamma) # estimate the BV frequency performing a center derivative

         N_square=zeros(length(P)-1)

@inbounds for i in 1:length(P)-1
              if (i==1)
        @fastmath N_square[i]=(g[i]+g[i+1])*(log(P[i+1]/d[i+1]^gamma) -log(P[i]/d[i]^gamma))/(2*gamma*(r[i+1]-r[i]))
              else
        @fastmath N_square[i]=g[i]*(log(P[i+1]/d[i+1]^gamma) - log(P[i-1]/d[i-1]^gamma))/(gamma*(r[i+1]-r[i-1]))
              end

           end

          return N_square
end
