function deadsimplegemv(TRANS::String, m::Int, n::Int, alpha::Float64, A::Array{MC{N},2}, x::Array{MC{N},1}, beta, y::Array{MC{N},1})
    if TRANS == "N"
        return alpha*A*x + beta*y
    else
        return alpha*transpose(A)*x + beta*y
    end
end
