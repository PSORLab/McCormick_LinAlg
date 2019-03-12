

function simplexscal(X::Array{MC{N},1}, a::Float64) where N
    R::Array{MC{N},1} = [a*x for x in X]
    return R
end

function deadsimplexscal(X::Array{MC{N},1}, a::Float64) where N
    R::Array{MC{N},1} = X*a
    return R
end
