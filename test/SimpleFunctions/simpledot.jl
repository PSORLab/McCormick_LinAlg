#Function to bechmark BLAS implementation against

function simpledot(X::SVector,Y::SVector)
    n = length(X)
    N = length(X[1].cv_grad)
    R::MC = MC{N}(0.,0.)
    for i in range(1,n)
        R += X[i] * Y[i]
    end
    return R
end

function deadsimpledot(X::SVector, Y::SVector)
    Z = SUM(X .* Y)
end
