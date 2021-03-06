 #replica of table 2 from DetofIntMatrices Paper (HORA´CEK)
 using IntervalArithmetic, LinearAlgebra
 include("SpecialFunctions.jl")
 include("CRAM.jl")
 include("CRAMinv.jl")
 include("GAUSS.jl")
 include("GAUSSinv.jl")
 include("GERSCH.jl")
 include("GERSCHinv.jl")
 include("HAD.jl")
 include("HADinv.jl")


function test()
 ns = [5,10,20,50]
 runs = 20
 lenn = length(ns)
 GE_inv::Array{Any, 1} = fill(0.0,lenn)
 GE_::Array{Any, 1} = fill(0.0,lenn)
 CRAM_inv::Array{Any, 1} = fill(0.0,lenn)
 GERSCH_inv::Array{Any, 1} = fill(0.0,lenn)

 GE_inv = fill(0.0,lenn)
 HAD_inv = fill(0.0,lenn)
 CRAM_inv = fill(0.0,lenn)
 GERSCH_inv = fill(0.0,lenn)


 lists = [GE_inv, HAD_inv, CRAM_inv, GERSCH_inv]
 labels = ["GE_inv", "HAD_inv", "CRAM_inv", "GERSCH_inv"] #GE not meant to be compared
 functs = [gaussinvdet, HADinv, craminvdet, GERSCHinv]
 m = length(lists)
 Inftol = 1E50
 r1 = .001
 r2 = .00001
 r1 = r2
for n = 1:lenn#For each dimension
      totals = fill(0, m) #Track runs for each method for average later
      for i = 1:runs #Perform all runs with dimesions ns[n] x ns[n]
            Am = Array{Float64,2}(undef,(ns[n],ns[n]))
            Am = map(x->rand()rand([-1,1]), Am)
            #=while det(Am) < 0  #Set conditions on matrices tested
                  Am = map(x->rand()rand([-1,1]), Am)
            end =#
            A_1 = map(x->Interval(x-r1, x+r1), Am)
            for j in 1:m #for each method tested
                  parts = functs[j](A_1)
                  #= if j == 1  #Check a specific estimation method
                        if diam(parts) > 1
                              println("here")
                              println(size(A_1))
                              println(parts)
                              println(det(Am))
                        end
                  end =#
                  part = diam(parts)
                  if part > Inftol || part < -Inftol || parts == ∅ || part == 0.0 #dont know why ∅ comes up. 0.0 case also confusing
                        #print(part)
                        continue
                  else
                        lists[j][n] += part
                        totals[j] += 1 #Don't count infinite ranges
                  end
            end
      end
      println(totals)
      #println(lists[:][n])
      for i = m:-1:1
            if totals[i] == 0
                  lists[i][n] = "No Runs" #Normalize by GEinv
            else
                  lists[i][n] = (lists[i][n] / totals[i]) /  (lists[1][n] / totals[1]) #Normalize by GEinv
            end
      end
end
return push!(lists, labels)
end
lists = test()
