

function simplesaxpy(a::Float64, X::SVector, Y::SVector)
    R::SVector{length(X)} = [a*x for x in X]
    R = [R[i] + Y[i] for i in range(1, length(X))]
    return R
end

function deadsimplesaxpy(a::Float64, X::SVector, Y::SVector)
    R::SVector{length(X)} = a*X .+ Y
    return R
end
