
using LinearAlgebra
function simplexnrm2(X::Array{MC{N},1}) where N
    R::MC = MC{N}(0.,0.)
    for x in X
        R += X*X
    end
    R = sqrt(R)
    return R
end

function deadsimplexnrm2(X::Array{MC{N},1}) where N
    Z::MC = LinearAlgebra.norm(X, 2)
    return Z
end
