

function simplexscal(a::Float64, X::SVector)
    R::SVector = [a*x for x in X]
end

function deadsimplexscal(a::Float64, X::SVector)
    R::SVector = X*a
end
