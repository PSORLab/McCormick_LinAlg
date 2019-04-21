

function GERSCH(A)
    Intzero::Interval = Interval(0.0,0.0)
    (n,m) = size(A)
    eigbounds = Array{Interval,1}(undef, n)
    for j = 1:n #Using Columns b/c of data storage
        middle = mid(A[j,j]) #Midpoint of disk
        rad = 0.0 #RAdius of disc
        for i = 1:n #Sum column to get disc radius
            if i == j
                continue
            end
            rad += mag(A[i,j])
        end
        eigbounds[j] = Interval(middle - rad, middle + rad)#1 is lower bounds interval, 2 is upper bound interval
    end
    detg::Interval = prod(eigbounds)
    return detg
end
