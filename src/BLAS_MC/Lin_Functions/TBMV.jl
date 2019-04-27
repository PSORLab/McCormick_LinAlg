#Triangular Banded Matrix vector product
#Not functional
function TBMV(UPLO::String, TRANS::String, DIAG::String, n::Integer, k::Integer, A::Array{MC{N},2}, LDA::Integer, x::Array{MC{N},1}) where N
#=
UPLO = 'U' or 'u'   Only the upper triangular part of A
     = 'L' or 'l'   Only the lower triangular part of A
DIAG = "U" if A is unit traingular (main diag is 1's)
     = "N" if not assumed
       =#
#Skip testing parameters
MCzero::MC = zero(MC{N})
#Not using sparse vectors, assume incx,incy==1
kx::Int =1
ky::Int =1
x_2::Array{MC{N},1} = x #Result vector. dont want to solve in place of right now
temp::MC = MCzero
nounit::Bool = (DIAG == "N")

if TRANS == "N" #A*x

    if UPLO == "U"
        kplus1 = k + 1
        for j in 1:n
            if x_2[j] != MCzero #May need to check this different way
                temp = x_2[j]
                l = kplus1 -j
                for i = max(1,j-k):(j-1)
                    x_2[i] += temp*A[l+i,j]
                end
                if nounit
                    x_2[j] *= A[kplus1,j]
                end
            end
        end

    else #Use lower triangular

        for j = n:-1:1
            if x_2[j] != MCzero
                temp = x_2[j]
                l = 1 - j
                for i = min(n,j+k):-1:(j+1)
                    x_2[j] += temp*A[l+i,j]
                end
                if nounit
                    x_2[j] *= A[1,j]
                end
            end
        end
    end

else #transpose(A)*x

    if UPLO == "U"
        kplus1 = k + 1
        for j = n:-1:1
            temp = x_2[j]
            l = kplus1 - j
            if nounit
                temp *= A[kplus1,j] #This line above may be i,j (wrong)
            end
            for i = (j-1):-1:max(1,j-k)
                temp += A[l+i,j]*x_2[i]
            end
            x_2[j] = temp
        end

    else

        for j = 1:n
            temp = x_2[j]
            l = 1 - j
            if nounit
                temp *= A[1,j]
            end
            for i = (j+1):min(n,j+k)
                temp += A[l+i,j]*x_2[i]
            end
            x_2[j] = temp
        end
    end
end
return x_2
end

#For SparseArrays
#Bascially a copy of GBMVs but not sure if it can get better
function TBMVs(TRANS::String, alpha::Float64, A::SparseMatrixCSC{MC{N},Int64}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N #A is a sparse matrix
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
    if TRANS == "N"
        for col = 1:A.n
            axj = alpha*x[col]
            for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                y_[rv[j]] += nzv[j]*axj
            end
        end
    else
        for col = 1:A.n
            temp = MCzero
            for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                temp += nzv[j]*x[rv[j]]
            end
            y_[col] += alpha * temp
        end

    end
    return y_
end
