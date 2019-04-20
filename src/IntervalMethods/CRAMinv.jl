
function craminvdet(A)#Lower bound seems too good to be true. Above HULL from paper
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
