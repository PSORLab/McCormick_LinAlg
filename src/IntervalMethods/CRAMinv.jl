
function solve(A_,b_) #Solve Ax=b
      (n,m) = size(A_)#Assumed to be square. Also no interval contains 0
      #println(n)
      A = copy(A_)
      b = copy(b_)
      #display(A)
      #display(b)
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

function CRAMinvdet(A)#Lower bound seems too good to be true. Above HULL from paper
Intzero::Interval = Interval(0.0,0.0)
(n,m) = size(A)
Am = map(x->(x.lo+x.hi)/2, A)# Midpoint Matrix for preconditioning
Am_inv = inv(Am)
detAm_inv = det(Am_inv)
A = A * Am_inv #Preconditioning step

x::Interval = Interval(1,1)#cumulative divisor
b::Array{Interval,1} = [Intzero]
Int_1 = Interval(1.,1.)
for d = 1:(n-1) #Number of dimensions left
      b = fill(Intzero, n-d+1)
      b[1] = Int_1
      x_d = solve(A[d:n, d:n], b)[1]
      x *= x_d
end
x_d = A[n,n] #Base case det(1x1 Matrix) = itself
detc = x_d / x
detc = detc / detAm_inv
return detc
end
