

function simplexscal(X::SVector, a::Float64)
    R::SVector{length(X)} = [a*x for x in X]
end

function deadsimplexscal(X::SVector, a::Float64)
    R::SVector{length(X)} = X*a
end
