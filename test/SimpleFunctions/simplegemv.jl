function deadsimplegemv(TRANS::String, m::Int, n::Int, alpha::Float64, A::Array{MC{N},2}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
    if TRANS == "N"
        return alpha*A*x + beta*y
    else
        return alpha*transpose(A)*x + beta*y
    end
end

#doesn't work. Cant mult A*x for unknown reason
function deadsimplegbmvs(TRANS::String, alpha::Float64, A::SparseMatrixCSC{MC{N},Int64}, x::Array{MC{N},1}, beta::Float64, y::Array{MC{N},1}) where N
    if TRANS == "N"
        return alpha*A*x + beta*y
    else
        return alpha*transpose(A)*x + beta*y
    end
end
