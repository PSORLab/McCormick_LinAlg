#Scalar multiplication of vector containing McCormick Objects
#MC{cc,cv,Intv,cc_grad,cv_grad,cnst}
#Implementation first found in Mitsos2009 (not full citation)
function XSCAL(MCv::SVector, scal::Float64)#where A<:AbstractArray #MC = single MC object, scal = Float64
    N::Int = length(MCv[1].cc_grad)
    n::Int = length(MCv)
    MCvtemp = Vector{MC}(undef, n)
    lo::Float64 = 0.0
    hi::Float64 = 0.0
    cc::Float64 = 0.0
    cc_grad = Vector{Float64}(undef, N)
    cv_grad = Vector{Float64}(undef, N)
    cnst::Bool = true
        if scal >= 0
                    for i in range(1,n)
                            x::MC = MCv[i]
                            lo = scal * x.Intv.lo
                            hi = scal * x.Intv.hi
                            cc = scal * x.cc
                            cv = scal * x.cv
                            cnst = x.cnst
                            for ip in range(1, N) #added length
                                        global cc_grad[ip] = scal * x.cc_grad[ip]
                                        global cv_grad[ip] = scal * x.cv_grad[ip]
                            end
                            global MCvtemp[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end
        else
                    for i in range(1,n)
                            x::MC = MCv[i]
                            lo = scal * x.Intv.hi
                            hi = scal * x.Intv.lo
                            cc = scal * x.cv
                            cv = scal * x.cc
                            cnst = x.cnst
                            for ip in range(1, N)
                                            global cc_grad[ip] = scal * x.cv_grad[ip]
                                            global cv_grad[ip] = scal * x.cc_grad[ip]
                            end
                            global MCvtemp[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end
        end
        return SVector{n, MC}(MCvtemp)
end
