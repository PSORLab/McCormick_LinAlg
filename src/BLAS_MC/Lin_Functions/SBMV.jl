#Symmetric Banded Matrix vector product
#Not functional
function SBMV(UPLO::String, n::Integer, k::Integer, alpha::Float64,  A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
 #ignoring parameter check
 #ignore sparse vectors so kx,ky = 1
 MCzero::MC = MC{N}(0.0, 0.0)
 #Not using sparse vectors, assume incx,incy==1
 kx::Integer =1
 ky::Integer =1
 #Form beta*y
y_2::Array{MC{N},1} = fill(MCzero, leny) #Result vector. dont want to solve in place of right now
 if beta != 1#If 1 can ignore coefficient
    if beta == 0
        y_2 .= MCzero
    else
        y_2[:] = XSCAL(beta, y)[:]#May replace this
    end
 end

 if alpha==0
     return y_2
 end

 temp1::MC = MCzero
 temp2::MC = MCzero
 l::Integer = 0

 if "U" == UPLO #use upper triangular of A
    kplus1 = k + 1
    for j = 1:n
         temp1 = alpha*x[j]
         temp2 = MCzero
         l = kplus1 - j
         for i = max(1,j-k):(j-1)
             y_2[i] += temp1*A[l+i,j]
             temp2 += A[l+i,j]*x[i]
         end
    y_2[j] += temp1*A[kplus1,j] + alpha*temp2
    end

else #Use lower tringular of A

    for j = 1:n
        temp1 = alpha*x[j]
        temp2 = MCzero
        y_2[j] += temp1*A[1,j]
        l = 1-j
        for i = (j+1):min(n,j+k)
            y_2[i] += temp1*A[l+i,j]
            temp2 += A[l+i,j]*x[i]
        end
        y_2[j] += alpha*temp2
    end
end
return y_2
end
