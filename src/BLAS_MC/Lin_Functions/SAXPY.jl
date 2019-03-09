#Scalar multiplication of vector containing McCormick Objects
#MC{fields of cv,cc,Intv,cv_grad,cc_grad,cnst}
function SAXPY(scal::Float64, X::SVector, Y::SVector)#where A<:AbstractArray #MC = single MC object, scal = Float64
    N::Int = length(X[1].cc_grad)
    n::Int = length(X)
    R = Vector{MC}(undef, n)
    lo::Float64 = 0.0
    hi::Float64 = 0.0
    cc::Float64 = 0.0
    cc_grad = Vector{Float64}(undef, N)
    cv_grad = Vector{Float64}(undef, N)
    cnst::Bool = true
        if scal >= 0
                    for i in range(1,n)
                            x = X[i]
                            y = Y[i]
                            lo = scal * x.Intv.lo + y.Intv.lo
                            hi = scal * x.Intv.hi + y.Intv.hi
                            cc = scal * x.cc      + y.cc
                            cv = scal * x.cv      + y.cv
                            cnst = x.cnst && y.cnst
                            for ip in range(1, N) #added length
                                        global cc_grad[ip] = scal * x.cc_grad[ip] + y.cc_grad[ip]
                                        global cv_grad[ip] = scal * x.cv_grad[ip] + y.cv_grad[ip]
                            end
                            global R[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end
        else
                    for i in range(1,n)
                            x = X[i]
                            y = Y[i]
                            lo = scal * x.Intv.hi + y.Intv.lo
                            hi = scal * x.Intv.lo + y.Intv.hi
                            cc = scal * x.cv      + y.cc
                            cv = scal * x.cc      + y.cv
                            cnst = x.cnst && y.cnst
                            for ip in range(1, N)
                                            global cc_grad[ip] = scal * x.cv_grad[ip] + y.cc_grad[ip]
                                            global cv_grad[ip] = scal * x.cc_grad[ip] + y.cv_grad[ip]
                            end
                            global R[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end
        end
        return SVector{n, MC}(R)
end
