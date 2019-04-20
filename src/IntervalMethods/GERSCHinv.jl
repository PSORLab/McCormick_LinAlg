
function GERSCHinv(A)
    Intzero::Interval = Interval(0.0,0.0)
    (n,m) = size(A)
    Am = map(x->(x.lo+x.hi)/2, A)# Midpoint Matrix for preconditioning
    Am_inv = inv(Am)
    detAm_inv = det(Am_inv)
    A = A * Am_inv #Preconditioning step
    circle = Array{Interval,1}(undef, 2)
    eigbounds = Array{Interval,1}(undef, n)
    for j = 1:n #Using Columns b/c of data storage
        circle = [A[j,j], A[j,j]]
        for i = 1:n
            if i == j
                continue
            end
            circle[:] = [circle[1]-abs(A[i,j]), circle[2]+ abs(A[i,j])]
        end
        eigbounds[j] = Interval(circle[1].lo, circle[2].hi)#1 is lower bounds interval, 2 is upper bound interval
    end
    detg::Interval = prod(eigbounds)
    detg = detg / detAm_inv
    return detg
end
