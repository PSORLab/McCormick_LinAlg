#Based on the BLAS implementation of GBMV

#Not including parameters LDA, INCX, INCY. So no sparse vectors.

function GBMV(TRANS::String, m::Integer, n::Integer, kl::Integer, ku::Integer, alpha::Float64,  A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
    #Ignoring dimensional check flags for now, assuming inputs are correct
    temp::MC = MC{N}(0.0, 0.0)
    if m==0 || n==0 || (alpha==0 && beta==0) #If returning real value, separate m,n==0 from alpha beta
        tempA::Array{MC{N},1} = [temp for ms in 1:m, ns in 1:n]
        return tempA
    end
    if TRANS == "N" #Do not use transpose
           lenx = n
           leny = m
    else#Use transpose
           lenx = m
           leny = n
    end
    y_2::Array{MC{N},1} = fill(MCzero, leny) #Result vector

    #form y_2 = beta*y
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
    kup1 = ku + 1
    #Vector Matrix multiply
    if TRANS == "N" #####Don't use transpose
        incx::Integer = 1
        jx::Integer = 1 #Assuming 1 == incx > 0
        for j in 1:n
            temp = alpha*x[jx]
            k = kup1 - j
            for i = max(1,j-ku):min(m,j+kl)
                y_2[i] = y_2[i] + temp*A[k+i,j]
            end
            jx = jx + incx
        end
    else############Use the transpose of A. Assume incx ==1
        incy::Integer = 1
        jy::Integer = 1 #Assuming 1 == incy > 0
        MCzero = temp
        for j in 1:n
            temp = MCzero
            k = kup1 - j
            for i in max(1, j-ku):min(m, j+kl)
                temp += A[k+i,j]*x[i]
            end
            y_2[jy] += alpha*temp
            jy += incy
        end
    end
    return y_2
end
