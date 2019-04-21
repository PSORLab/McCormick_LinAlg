#Not perfect. Also useless because it does not guarentee a sign for the det.
function HADinv(A::Array{Interval{Float64},2})
    Am = map(x->(x.lo+x.hi)/2, A)# Midpoint Matrix for preconditioning
    Am_inv = inv(Am)
    detAm_inv = det(Am_inv)
    A = Am_inv * A #Preconditioning step

    (m,n) = size(A)
    Int_1 = Interval(1.,1.)
    deth::Interval = Interval(1.,1.)
    colnorm::Interval = Int_1
    for j = 1:n
        colnorm = Int_1
        for i = 1:n #Assume square matrix
            colnorm += A[i,j]^2
        end
        colnorm = sqrt(colnorm) #Is this a defined function for IntervalArith?
        deth = deth * colnorm
    end
    d = max(abs(deth.lo), abs(deth.hi))#Take max abs value of interval for bounds on deth
    deth = Interval(-d, d) / detAm_inv
    return deth
end
