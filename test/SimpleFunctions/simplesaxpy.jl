

function simplesaxpy(a::Float64, X::Array{MC{N},1}, Y::Array{MC{N},1}) where N
    R::Array{MC{N},1} = [a*x for x in X]
    R = [R[i] + Y[i] for i in range(1, length(X))]
    return R
end

function deadsimplesaxpy(a::Float64, X::Array{MC{N},1}, Y::Array{MC{N},1}) where N
    R::Array{MC{N},1} = a*X .+ Y
    return R
end
