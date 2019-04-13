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
function GEMV(TRANS::String, m::Int, n::Int, alpha::Float64, A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
#test of parameter correctness skipped\
    temp::MC = MC{N}(0.0,0.0)#Could write this out so not to outsource init
    MCzero::MC = MC{N}(0.0,0.0)
    if TRANS == "N" #Do not use transpose
           lenx = n
           leny = m
    else#Use transpose
           lenx = m
           leny = n
    end
    incx::Int, incy::Int = 1,1
    y_2::Array{MC{N},1} = fill(MCzero, leny) #Result vector
    ###############form y_2 = beta*y
    if beta != 1#If 1 can ignore coefficient
    if beta == 0
        y_2 .= temp
    else
        y_2[:] = XSCAL(beta, y)[:]#May replace this
    end
    end

    if alpha == 0.0 # Second term disappears
        return y
    end


    if TRANS == "N"
        #############make y := alpha*A*x + y
        jx::Int = 1
        for j in 1:n
            temp = alpha*x[jx]
            for i in 1:m
                y_2[i] = y_2[i] + temp*A[i,j]
            end
            jx = jx + incx
            #if incx: jx = jx + incx used for vector temp index
        end
    else #for y := alpha*A**T*x + y
        jy::Int = 1
        for j in 1:n
            temp = MCzero
            for i in 1:m
                temp += A[i,j] * x[i]
            end
            y_2[jy] = y_2[jy] + alpha*temp
            jy+= incy
            #With incy: inc jy by incy
        end
    end
return y_2
end
