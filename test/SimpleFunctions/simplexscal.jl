

function simplexscal(a::Float64, X::Array{MC{N},1}) where N
    R::Array{MC{N},1} = [a*x for x in X]
    return R
end

function deadsimplexscal(a::Float64, X::Array{MC{N},1}) where N
    R::Array{MC{N},1} = X*a
    return R
end
