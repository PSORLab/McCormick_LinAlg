using IntervalArithmetic
#Want to make sure type correct functions are called, not generic ones

function dual(Inter::Interval) # E. Kaucher, Interval analysis in the extended Interval space IR, Computing, Suppl., 2 (1980), 33-49, doi: 10.1007/978-3-7091-8577-3
      return Interval(Inter.hi, Inter.lo)
end
function dualsub(Inter1, Inter2)
      #print("Generic Method")
      return Interval(Inter1.lo - Inter2.lo, Inter1.hi-Inter2.hi)
end
function dualsub(Inter1::Interval, Inter2::Interval)
      return Interval(Inter1.lo - Inter2.lo, Inter1.hi-Inter2.hi)
end
function dualsub(Inter1::Array{Interval{Float64},1}, Inter2::Interval)
      return map(x->dualsub(x, Inter2), Inter1)
end
function dualsub(Inter1::Array{Interval{Float64},1}, Inter2::Array{Interval{Float64},1})
      return map(dualsub, Inter1, Inter2)
end

function dualdiv(Inter1, Inter2)
      #print("Generic Method")
      return Inter1 / dual(Inter2)
end
function dualdiv(Inter1::Interval, Inter2::Interval)
      return Inter1 / dual(Inter2)
end
function dualdiv(Inter1::Array{Interval,1}, Inter2::Interval)
      return Inter1 / dual(Inter2)
end
