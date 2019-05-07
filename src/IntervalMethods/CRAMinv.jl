
function craminvdet(A::Array{Interval{Float64}, 2})#Lower bound seems too good to be true. Above HULL from paper
Intzero::Interval = Interval(0.0,0.0)
(n,m) = size(A)
Am = map(x->(x.lo+x.hi)/2, A)# Midpoint Matrix for preconditioning
Am_inv = inv(Am)
detAm_inv = det(Am_inv)
A = Am_inv * A  #Preconditioning step

x::Interval = Interval(1,1)#cumulative divisor
b::Array{Interval{Float64},1} = [Intzero]
Int_1 = Interval(1.,1.)
for d = 1:(n-1) #Number of dimensions left
      b = fill(Intzero, n-d+1)
      b[1] = Int_1
      x_d = solve2(A[d:n, d:n], b)[1]
      x *= x_d
end
x_d = A[n,n] #Base case det(1x1 Matrix) = itself
detc = x_d / x
detc = detc / detAm_inv
return detc
end

function craminvdet(A::Array{MC{N},2}) where N#Lower bound seems too good to be true. Above HULL from paper
Intzero::Interval = Interval(0.0,0.0)
(n,m) = size(A)
B = map(x -> x.Intv, A)
Am = map(x->(x.lo+x.hi)/2, B)# Midpoint Matrix for preconditioning
Am_inv = inv(Am)
detAm_inv = det(Am_inv)
A = Am_inv * B  #Preconditioning step

x::Interval = Interval(1,1)#cumulative divisor
b::Array{Interval{Float64},1} = [Intzero]
Int_1 = Interval(1.,1.)
for d = 1:(n-1) #Number of dimensions left
      b = fill(Intzero, n-d+1)
      b[1] = Int_1
      x_d = solve2(B[d:n, d:n], b)[1]
      x *= x_d
end
x_d = B[n,n] #Base case det(1x1 Matrix) = itself
detc = x_d / x
detc = detc / detAm_inv
return detc
end
