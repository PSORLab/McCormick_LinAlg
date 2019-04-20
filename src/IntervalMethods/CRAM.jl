
function cramdet(A::Array{Interval{Float64},2})#Lower bound seems too good to be true. Above HULL from paper
Intzero::Interval = Interval(0.0,0.0)
(n,m) = size(A)

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
return detc
end

function cramdet(A::Array{Float64,2})#Lower bound seems too good to be true. Above HULL from paper
(n,m) = size(A)
x::Float64 = 1.
b::Array{Float64,1} = [0.]
for d = 1:(n-1) #Number of dimensions left
      b = fill(0, n-d+1)
      b[1] = 1
      x_d = solve(A[d:n, d:n], b)[1] #This solve is wrong
      #x_d = (inv(A[d:n,d:n])* b)[1]
      x *= x_d
end
x_d = A[n,n] #Base case det(1x1 Matrix) = itself
detc = x_d / x
return detc
end
