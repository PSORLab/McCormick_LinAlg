#Function to bechmark BLAS implementation against

function simpledot(X::Array{MC{N},1},Y::Array{MC{N},1}) where N
    n = length(X)
    R::MC = MC{N}(0.,0.)
    for i in range(1,n)
        R += X[i] * Y[i]
    end
    return R
end

function deadsimpledot(X::Array{MC{N},1}, Y::Array{MC{N},1}) where N
    Z::MC = sum(X .* Y)
    return Z
end
