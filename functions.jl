# -------------------------------------------------------------------------
function compute_price_tradeshares(w::Array{Float64,1},
                                   A::Array{Float64,1},
                                   d::Array{Float64,2},
                                   theta::Float64,
                                   gam::Float64,
                                   N::Int64,
                                   v::Float64,
                                   P0::Array{Float64,1})

   #  Predifne some variables
   bts = Array{Float64, 2}(undef, N,N)
   P = Array{Float64, 1}(undef, N)

   P = copy(P0)
   Tp = copy(P0)
   u = Array{Float64, 1}(undef, N) #???

   tol = 1.0e-3 # Tolerance for convergence
   dif = 10.0 # Initialize difference from prices
   iter = 0

   while (dif>tol) & (iter<10000)
       iter += 1

       #P = copy(Tp)

       for n in 1:N
           u[n] = ((w[n]./v).^(v) .* (P[n]./(1.0.-v)).^(1.0.-v)) #????
       end #For n
       #update tp
       for n in 1:N
           denom = 0.0
           for i in 1:N
               denom += (d[n,i].^(-theta) .* ((u[i].^(-theta)).*(A[i]).^(v.*theta)))
           end # for2

           Tp[n] = denom.^(-1.0./theta)

           for i in 1:N
               bts[n,i] = d[n,i].^(-theta) .* ((u[i].^(-theta)).*(A[i]).^(v.*theta)) ./
                          denom
           end #for
       end

       #pmax = (Tp .- P)./(Tp)
       pmax = Tp .- P
       dif = maximum(abs.(pmax)) #?????
       scl = 0.2
       P = P .+ scl.* pmax

      @printf " %d iterations. price diff: %f \n" iter dif

       # Compute price vector and trade share matrix

      end #while

    return P, bts, u
end
# -------------------------------------------------------------------------




# -------------------------------------------------------------------------
function compute_wage(w0::Array{Float64,1},
                      A::Array{Float64,1},
                      d::Array{Float64,2},
                      theta::Float64,
                      gam::Float64,
                      N::Int64,
                      v::Float64,
                      P0::Array{Float64,1})



  #  Predifne some variables
  bts = Array{Float64, 2}(undef, N,N)
  P = Array{Float64, 1}(undef, N)
  EX = Array{Float64, 1}(undef, N)
  IM = Array{Float64, 1}(undef, N)
  u = Array{Float64, 1}(undef, N)
  w = copy(w0)
  Tw = copy(w0)


  tol = 1.0e-4 # Tolerance for convergence
  dif = 10.0 # Initialize difference from trade balance
  iter = 0
  while (dif>tol) & (iter<10000)
      iter += 1

      w = copy(Tw)

      P, bts, u = compute_price_tradeshares(w,
                                         A,
                                         d,
                                         theta,
                                         gam,
                                         N,
                                         v,
                                         P0)

    
    #   tot_prod = Array{Float64,1}(undef, N);
    #   tot_prod = w.*L./v;

    #   bts1 = transpose(bts);
    #   tot_exp = Array{Float64, 1}(undef, N);
    #   tot_exp = bts1\tot_prod;

    #    P0 = copy(P)
      # Compute each country's exports and imports
      EX = Array{Float64, 1}(undef, N)
      IM = Array{Float64, 1}(undef, N)
      for n in 1:N
          EX[n] = 0
          IM[n] = 0
          for i in 1:N
              if i != n
                  EX[n] += w[i].*L[i].*bts[i,n]
                  IM[n] += w[n].*L[n].*bts[n,i]

                #   IM[n] += tot_exp[n].*L[n].*bts[n,i]
                #   EX[n] += tot_exp[i].*L[i].*bts[i,n]
              end # end if
          end # for i
      end # for n

      # Define excess demand (net exports to GDP)
      Zw = (EX .- IM)./(w.*L)
      dif = maximum(abs.(Zw))

      # Update the guess for wages
      scl = 0.1
      Tw = w.*(1.0 .+ scl.*Zw) # new, updated guess for w

      if mod(iter,100)==0
         # @printf " %d iterations. Excess demand: %f \n" iter dif
      end

  end # while

  return P, bts, u, w, EX, IM

end
# -------------------------------------------------------------------------
