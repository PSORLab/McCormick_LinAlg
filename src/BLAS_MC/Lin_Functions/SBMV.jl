#Symmetric Banded Matrix vector product
function SBMV(UPLO::String, n::Integer, k::Integer, alpha::Float64,  A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
 #ignoring parameter check
 #ignore sparse vectors so kx,ky = 1
 MCzero::MC = zero(MC{N})
 #Not using sparse vectors, assume incx,incy==1
 kx::Int =1
 ky::Int =1
 #Form beta*y
y_2::Array{MC{N},1} = Array{MC{N},1}(undef, length(y)) #Result vector. dont want to solve in place of right now
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
 l::Int = 0

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

#Still a copy of GBMV right now. Can it get better though if in full storage?
function SBMVs(alpha::Float64, A::SparseMatrixCSC{MC{N},Int64}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N #A is a sparse matrix
    nzv = A.nzval
    rv = A.rowval
    if TRANS == "N" #Do not use transpose
           lenx = A.n
           leny = A.m
    else#Use transpose
           lenx = A.m
           leny = A.n
    end
    y_::Array{MC{N},1} = copy(y)
    if beta != 1
        if beta != 0
            y_ = XSCAL!(beta, y_)
        else
            y_ = fill(y_, zero(eltype(y_)))
        end
    end
    #println(y_)
#    if TRANS == "N"
        for col = 1:A.n
            axj = alpha*x[col]
            for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                y_[rv[j]] += nzv[j]*axj
            end
        end
#    else
#=        for col = 1:A.n
            temp = MCzero
            for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                temp += nzv[j]*x[rv[j]]
            end
            y_[col] += alpha * temp
        end
        =#
#    end
    return y_
end
