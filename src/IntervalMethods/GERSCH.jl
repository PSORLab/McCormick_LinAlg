

function GERSCH(A)
    Intzero::Interval = Interval(0.0,0.0)
    (n,m) = size(A)
    circle = Array{Interval,1}(undef, 2)
    eigbounds = Array{Interval,1}(undef, n)
    for j = 1:n #Using Columns b/c of data storage
        circle = [A[j,j], A[j,j]]
        for i = 1:n #Sum column to get disc radius
            if i == j
                continue
            end
            circle[:] = [circle[1]-abs(A[i,j]), circle[2]+ abs(A[i,j])] #Each bound is center +/- radius
        end
        eigbounds[j] = Interval(circle[1].lo, circle[2].hi)#1 is lower bounds interval, 2 is upper bound interval
    end
    detg::Interval = prod(eigbounds)
    return detg
end
