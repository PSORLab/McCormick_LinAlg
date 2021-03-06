#Symmetric Matrix vector product
#Not functional
function SYMV(UPLO::String, n::Integer, alpha::Float64,  A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
#=
  UPLO is CHARACTER*1
    UPLO specifies whether the upper or lower
    triangular part of the array A is to be referenced as
    follows:
       UPLO = 'U' or 'u'   Only the upper triangular part of A
       is to be referenced.
       UPLO = 'L' or 'l'   Only the lower triangular part of A
       is to be referenced.
       =#
#Skip testing parameters
MCzero::MC = zero(MC{N})
#Not using sparse vectors, assume incx,incy==1
kx::Int =1
ky::Int =1
temp1::MC = MCzero
#Form beta*y
y_2::Array{MC{N},1} = Array{MC{N},1}(undef, length(y)) #Result vector. dont want to solve in place of right now
if beta != 1#If 1 can ignore coefficient
if beta == 0
    y_2 .= temp1
else
    y_2[:] = XSCAL(beta, y)[:]#May replace this
end
end

if alpha==0
    return y_2
end
temp2::MC = MCzero

if "U" == UPLO #Use upper triangular of A
    for j = 1:n
        temp1 = alpha*x[j]
        temp2 = MCzero
        for i = 1:(j-1)
            y_2[i] += temp1*A[i,j]
            temp2 += A[i,j]*x[i]
        end
        y_2[j] += temp1*A[j,j] + alpha*temp2
    end
else #Use lower triangular of A
    for j = 1:n
        temp1 = alpha*x[j]
        temp2 = MCzero
        y_2[j] += temp1*A[j,j]
        for i = (j+1):n
            y_2[i] += temp1*A[i,j]
            temp2 += A[i,j]*x[i]
        end
        y_2[j] += alpha*temp2
    end
end
return y_2
end
