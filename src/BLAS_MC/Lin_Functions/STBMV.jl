#Triangular Banded Matrix vector product
#Not functional
function STBMV(UPLO::String, TRANS::String, DIAG::STRING, n::Integer, k::Integer, A::Array{MC{N},2}, LDA::Integer, x::Array{MC{N},1}) where N
#=
UPLO = 'U' or 'u'   Only the upper triangular part of A
     = 'L' or 'l'   Only the lower triangular part of A
DIAG = "U" if A is unit traingular (main diag is 1's)
     = "N" if not assumed
       =#
#Skip testing parameters
MCzero::MC = MC{N}(0.0, 0.0)
#Not using sparse vectors, assume incx,incy==1
kx =1
ky =1
x_2::Array{MC{N},1} = x #Result vector. dont want to solve in place of right now
temp::MC = MCzero
nounit = (DIAG == "N")

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
