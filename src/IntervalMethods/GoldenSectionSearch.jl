#Originally from an open source Python implementation for GSS found on https://en.wikipedia.org/wiki/Golden-section_search

#Does not reuse function evaluations, returns value
function gss(f, a, b, tol=1e-5)
    gr = (sqrt(5) + 1) / 2
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(c - d) > tol
        if f(c) < f(d)
            b = d
        else
            a = c
        end
        # we recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
    end

    return (b + a) / 2
end


#Reuses Function evaluations, returns interval
function gss2(f,a,b,tol=1e-5)
    invphi = (sqrt(5) - 1) / 2 # 1/phi
    invphi2 = (3 - sqrt(5)) / 2 # 1/phi^2

    (a,b)=(min(a,b),max(a,b))
    h = b - a
    if h <= tol
         return (a,b)
    end
    # required steps to achieve tolerance
    n = Int(ceil(log(tol/h)/log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for k in 1:(n-1)
        if yc < yd
            b = d
            d = c
            yd = yc
            h = invphi*h
            c = a + invphi2 * h
            yc = f(c)
        else
            a = c
            c = d
            yc = yd
            h = invphi*h
            d = a + invphi * h
            yd = f(d)
        end
    end

    if yc < yd
        return (a,d)
    else
        return (c,b)
    end
end
