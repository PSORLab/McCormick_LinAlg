#Scalar multiplication of vector containing McCormick Objects
#MC{cc,cv,Intv,cc_grad,cv_grad,cnst}
#Implementation first found in Mitsos2009 (not full citation)
function XSCAL(MCv::SVector{n, MC{N}}, scal::Float64) where n <: Integer where N <: Integer#where A<:AbstractArray #MC = single MC object, scal = Float64
                if scal >= 0
                                for x in MCv
                                                x.lo = scal * x.lo
                                                x.hi = scal * x.hi
                                                x.cc = scal * x.cc
                                                x.cv = scal * x.cv
                                                for ip in range(1, length(MC.cc_grad)) #added length
                                                                x.cc_grad[ip] = scal * x.cc_grad[ip]
                                                                x.cv_grad[ip] = scal * x.cv_grad[ip]
                                                end
                                end
                else
                                for x in MCv
                                                x.lo = scal * x.hi
                                                x.hi = scal * x.lo
                                                x.cc = scal * x.cv
                                                x.cv = scal * x.cc
                                                for ip in range(1, length(MC.cc_grad))
                                                                MC.cc_grad[ip] = scal * MC.cv_grad[ip]
                                                                MC.cv_grad[ip] = scal * MC.cc_grad[ip]
                                                end
                                end
                end
                return MCv
end
