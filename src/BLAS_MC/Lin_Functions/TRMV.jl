#Triangular Matrix vector product
#Not functional
function TRMV(UPLO::String, TRANS::String, DIAG::String, n::Integer, A::Array{MC{N},2}, x::Array{MC{N},1}) where N
#=
UPLO = 'U' or 'u'   Only the upper triangular part of A
     = 'L' or 'l'   Only the lower triangular part of A
DIAG = "U" if A is unit traingular (main diag is 1's)
     = "N" if not assumed
       =#
#Skip testing parameters
MCzero::MC = MC{N}(0.0, 0.0)
#Not using sparse vectors, assume incx,incy==1
kx::Int =1
ky::Int =1
leny = length(x)
x_2::Array{MC{N},1} = fill(MCzero, leny)  #Result vector
for i in 1:leny
    x_2[i] = x[i]
end
temp::MC = MCzero
nounit = (DIAG == "N")
if TRANS == "N" #A*x

    if UPLO == "U"#Upper triangular

        for j in 1:n
            if x_2[j] != MCzero #May need to check this different way
                temp = x_2[j]
                for i = 1:(j-1)
                    x_2[i] += temp*A[i,j]
                end
                if nounit
                    x_2[j] *= A[j,j]
                end
            end
        end

    else

        for j = n:-1:1 #Lower triangular
            if x_2[j] != MCzero
                temp = x_2[j]
                for i = n:-1:(j+1)
                    x_2[i] += temp*A[i,j]
                end
                if nounit
                    x_2[j] *= A[j,j]
                end
            end
        end
    end

else #For the transpose of A case

    if UPLO == "U" #Upper triangular

        for j = n:-1:1
            temp = x_2[j]
            if nounit
                temp *= A[j,j] #This line above may be i,j (wrong)
            end
            for i = (j-1):-1:1
                temp += A[i,j]*x_2[i]
            end
            x_2[j] = temp
        end

    else

        for j = 1:n #Lower triangular
            temp = x_2[j]
            if nounit
                temp *= A[j,j]
            end
            for i = (j+1):n
                temp += A[i,j]*x_2[i]
            end
            x_2[j] = temp
        end
    end
end
return x_2
end
