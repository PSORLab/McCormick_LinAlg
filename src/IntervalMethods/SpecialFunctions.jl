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
#Linear Solver using forward+back substitution
function solve(A_::Array{Interval{Float64},2},b_::Array{Interval{Float64},1}) #Solve Ax=b
      (n,m) = size(A_)
      A = copy(A_)
      b = copy(b_)
      Intzero::Interval = Interval(0.0,0.0)
      x::Array{Interval, 1} = fill(Intzero, n)
      for i = 1:n#Forward sub
            for j = (i+1):n
                  m = A[j,i]/A[i,i]
                  for k = i:n
                        A[j,k] = A[j,k] - m * A[i,k]
                        b[j] = b[j] - m * b[i]
                  end
            end
      end
      s::Interval = Intzero
      for i = 1:n #Backwards sub
            s = Intzero
            for j = 1:n
                  s += A[i,j]*x[j]
            end
            x[i] = (b[i]-s) / A[i,i]
      end

      return x
end
function solve(A_::Array{Float64,2}, b_::Array{Float64,1}) #Solve Ax=b
      (n,m) = size(A_)
      A = copy(A_)
      b = copy(b_)
      x::Array{Float64, 1} = fill(0, n)
      for i = 1:n#Forward sub
            for j = (i+1):n
                  m = A[j,i]/A[i,i]
                  for k = i:n
                        A[j,k] = A[j,k] - m * A[i,k]
                        b[j] = b[j] - m * b[i]
                  end
            end
      end
      s::Float64 = 0.0
      for i = 1:n #Backwards sub
            s = 0.0
            for j = 1:n
                  s += A[i,j]*x[j]
            end
            x[i] = (b[i]-s) / A[i,i]
      end

      return x
end

function solve2(A,b)
      A_ = copy(A)
      b_ = copy(b)
      Am = map(x->(x.lo+x.hi)/2, A_)
      bm = map(x->(x.lo+x.hi)/2, b_)

      y = inv(Am) * bm #Approx solution
      Az = b_ - A_*y #Solve for z, may need part 1 of this paper
end
