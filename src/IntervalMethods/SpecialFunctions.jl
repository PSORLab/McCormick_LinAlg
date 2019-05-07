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
      #Used by CRAMinv. Should make version that makes lower triangular to see if bounds on x[1] is different
function solve(A_::Array{Interval{T},2},b_::Array{Interval{T},1}) where T<:Real #Solve Ax=b
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
                  end
                  b[j] = b[j] - m * b[i]
            end
      end
      s::Interval = Intzero
      for i = n:-1:1 #Backwards sub
            s = Intzero
            for j = 1:n
                  s += A[i,j]*x[j]
            end
            #=println("Here")
            display(x[i])
            display(b[i])
            display(A[i,i]) =#
            x[i] = (b[i]-s) / A[i,i]
      end

      return x
end
function solve2(A_::Array{Interval{T},2},b_::Array{Interval{T},1}) where T<:Real #Solve Ax=b
      (m, n) = size(A_);
      A_kplus1 = copy(A_)
      A_k = copy(A_)
      b_kplus1 = copy(b_)
      b_k = copy(b_)
      Intzero = Interval(0.0, 0.0)
      for k = 1:(n-1) #From A NEW CRITERION TO GUARANTEE THE FEASIBILITY OF THE INTERVAL GAUSSIAN ALGORITHM*, A. FROMMERf AND G. MAYER
            for i = 1:n #Consider Pivot Tightening?
                  for j = 1:n#Assume A is nxn
                        if 1<=i<=k && 1<=j<=n
                        elseif (k+1)<=i<=n && (k+1)<=j<=n
                              A_kplus1[i,j] = A_k[i,j] - (A_k[i,k]*A_k[k,j])/A_k[k,k]
                              b_kplus1[i] = b_k[i] - (A_k[i,k]*b_k[k])/A_k[k,k]
                        else
                              A_kplus1[i,j] = Intzero
                        end
                  end
            end
            A_k = A_kplus1
            b_k = b_kplus1
      end
      x::Array{Interval{T}, 1} = copy(b_k)
      for i = n:-1:1
            s = Intzero
            for j = (i+1):n
                  s += A_k[i,j]*x[j]
            end
            x[i] = (b_k[i] - s)/A_k[i,i]
      end
      return x
end

function solve(A_::Array{Float64,2}, b_::Array{Float64,1}) #Solve Ax=b
      (n,m) = size(A_)
      A = copy(A_)
      b = copy(b_)
      display(b)
      x::Array{Float64, 1} = fill(0, n)
      for i = 1:n#Forward sub
            for j = (i+1):n
                  m = A[j,i]/A[i,i]
                  for k = i:n
                        A[j,k] = A[j,k] - m * A[i,k]
                  end
                  b[j] = b[j] - m * b[i]
            end
      end
      s::Float64 = 0.0
      for i = n:-1:1 #Backwards sub
            s = 0.0
            for j = 1:n
                  s += A[i,j]*x[j]
            end
            x[i] = (b[i]-s) / A[i,i]
      end

      return x
end
# m(a) = mid(a)), w = width = diam(a).
function intadd(a::Interval{T}, b::Interval{T}) where T<: Real #Real numbers should be treated as degenerate interval
      mab = mid(a) + mid(b)
      k = ((b.hi+a.hi) - (b.lo+a.lo))/2
      return Interval(mab-k, mab + k)
end
function intsub(a::Interval{T}, b::Interval{T}) where T<: Real
      if a == b
            #a-b = a-dual(a) = [0.0]
            return Interval(0.0,0.0)
      else
            mab = mid(a) - mid(b)
            k = ((b.hi+a.hi) - (b.lo+a.lo))/2
            return Interval(mab-k, mab + k)
      end
end
function intsub(Inter1::Array{Interval{T},1}, Inter2::Array{Interval{T},1}) where T <: Real
      result::Array{Interval{T},1} = []
      for i in 1:length(Inter1)
            push!(result, intsub(Inter1[i], Inter2[i]))
      end
      return result
end
function intmult(a::Interval{T}, b::Interval{T}) where T<: Real
      ab = a*b
      mab = mid(a) * mid(b)
      (alpha, beta) = (ab.lo, ab.hi)
      k = min(mab - alpha, beta - mab)
      return Interval(mab-k, mab + k)
end
function intmult(Inter1::Array{Interval{T},1}, Inter2::Interval{T}) where T <: Real
      result::Array{Interval{T},1} = []
      for i in Inter1
            push!(result, intmult(i, Inter2))
      end
      return result
end
function intdiv(a::Interval{T}, b::Interval{T}) where T<: Real
      if b.lo*b.hi <= 0
            return "Zero In Denominator"#Should raise error
      end
      if a == b
            #a/b = a/dual(a) = 1
            return Interval(1,1)
      else
            blo, bhi = b.lo, b.hi
            k = min((1/bhi)*((bhi-blo)/(blo+bhi)), (1/blo)*((bhi-blo)/(blo+bhi)))
            return intmult(a, Interval((1/mid(b))-k, (1/mid(b)) + k))
      end
end
function intdiv(Inter1::Array{Interval{T},1}, Inter2::Interval{T}) where T <: Real
      result::Array{Interval{T},1} = []
      for i in Inter1
            push!(result, intdiv(i, Inter2))
      end
      return result
end
