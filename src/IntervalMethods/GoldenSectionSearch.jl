#This is an open source implementation for GSS found on https://en.wikipedia.org/wiki/Golden-section_search

#=Python program for golden section search.  This implementation
   reuses function evaluations, saving 1/2 of the evaluations per
   iteration, and returns a bounding interval.=#
invphi = (math.sqrt(5) - 1) / 2 # 1/phi
invphi2 = (3 - math.sqrt(5)) / 2 # 1/phi^2
def gss(f,a,b,tol=1e-5):
    #=
    Given a function f with a single local minimum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the minimum with d-c <= tol.
    =#

    (a,b)=(min(a,b),max(a,b))
    h = b - a
    if h <= tol: return (a,b)

    # required steps to achieve tolerance
    n = int(math.ceil(math.log(tol/h)/math.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for k in xrange(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi*h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi*h
            d = a + invphi * h
            yd = f(d)

    if yc < yd:
        return (a,d)
    else:
        return (c,b)
