#Not perfect. Also useless because it does not guarentee a sign for the det.
function HADinv(A::Array{Interval{Float64},2})
    Am = map(x->(x.lo+x.hi)/2, A)# Midpoint Matrix for preconditioning
    Am_inv = inv(Am)
    detAm_inv = det(Am_inv)
    A = Am_inv * A #Preconditioning step
    (m,n) = size(A)
    Int_1 = Interval(1.,1.)
    colnorm::Float64 = 1.0
    d::Float64 = 1.0
    for j = 1:n #Product of Interval Eucl Norm of each column
        colnorm = 1.0
        for i = 1:n #Assume square matrix
            colnorm += mag(A[i,j])^2
        end
        colnorm = sqrt(colnorm)
        d *= colnorm
    end
    deth::Interval = Interval(-d,d) / detAm_inv
    return deth
end
