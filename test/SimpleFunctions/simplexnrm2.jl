
using LinearAlgebra
function simplexnrm2(X::SVector)
    R::MC = MC{length(X[1].cv_grad)}(0.,0.)
    for x in X
        R += X*X
    end
    R = sqrt(R)
    return R
end

function deadsimplexnrm2(X::SVector)
    LinearAlgebra.norm(X, 2)
end
