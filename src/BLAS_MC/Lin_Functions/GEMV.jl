#Based on the BLAS implementation of GEMV
#Performs alpha*A*x + beta*y if TRANS=N OR alpha*A**Tx + beta*y if TRANS=T
#=PARAMETERS
alpha, beta ::Float64, scalar coeff of a matrix
M,N ::Int, Matrix Dimensions
A ::Array, m x n Matrix
x,y ::Vectors
Not Implemented:
INCX, INCY ::Int!=0 for increment of x+y
LDA ::Int
=#
#Might actually set!(y) param to result
#Still need to write out field additions
function GEMV(TRANS::String, m::Int, n::Int, alpha::Float64, A::Array{MC{N},2}, x::Array{MC{N},1}, beta, y::Array{MC{N},1}) where N
#test of parameters skipped
    if TRANS == "N"
           lenx = n
           leny = m
    else
           lenx = m
           leny = n
    end
#no incx incy terms
    #form y = beta*y
    if beta != 1
    if beta == 0
        MCzero = MC{}
        y.= MCzero
    else
        y[:] = XSCAL(beta, y)[:]#May replace this
    end
    end
    if alpha == 0.0
        return y
    end
    temp::MC = MC{1}(0.0,0.0)
    if TRANS == "N"
        #For y := alpha*A*x + y
        for j in 1:n
            temp = alpha*x[j]
            for i in 1:m
                y[i] = y[i] + temp*A[i,j]
            end
            #if incx: jx = jx + incx used for vector temp index
        end
    else #for y := alpha*A**T*x + y
        for j in 1:n
            temp = MC{N}(0.0,0.0)
            for i in 1:m
                temp += A[i,j] * x[i]
            end
            y[j] = y[j] + alpha*temp
            #With incy: inc jy by incy
        end
    end
return y
end
