#Scalar multiplication of vector containing McCormick Objects
#MC{fields of cv,cc,Intv,cv_grad,cc_grad,cnst}
#Still needs mod5 added
function SAXPY(scal::Float64, X::Array{MC{N},1}, Y::Array{MC{N},1}) where N #where A<:AbstractArray #MC = single MC object, scal = Float64
    #N::Int = length(X[1].cc_grad)
    n::Int = length(X)
    R = Vector{MC}(undef, n) #Still needs the mod5 implementation
    lo::Float64 = 0.0
    hi::Float64 = 0.0
    cc::Float64 = 0.0
    cc_grad = Vector{Float64}(undef, N)
    cv_grad = Vector{Float64}(undef, N)
    cnst::Bool = true
    x::MC = MC{1}(0.0,0.0)
    y::MC = copy(x)

        if scal >= 0
                    for i in 1:n #using eachindex(x) or something else to create iterable? 1:n seems faster than range
                            x = X[i] #referecne list everytime, inbounds macro, first loop outside to preallo
                            y = Y[i]#look to julia performance tips, fastmath has slightly higher rounding error.@warrant_type
                            lo = scal * x.Intv.lo + y.Intv.lo #Just take and pass regular vectors. PS static ~<100-1000
                            hi = scal * x.Intv.hi + y.Intv.hi
                            cc = scal * x.cc      + y.cc
                            cv = scal * x.cv      + y.cv
                            cnst = x.cnst && y.cnst
                            for ip in 1:N #added length. Do same 5peat here?
                                        cc_grad[ip] = scal * x.cc_grad[ip] + y.cc_grad[ip]
                                        cv_grad[ip] = scal * x.cv_grad[ip] + y.cv_grad[ip]
                            end
                            R[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end #global not necessary in f(x)s!!!
        else
                    for i = 1:n
                            x = X[i]
                            y = Y[i]
                            lo = scal * x.Intv.hi + y.Intv.lo
                            hi = scal * x.Intv.lo + y.Intv.hi
                            cc = scal * x.cv      + y.cc
                            cv = scal * x.cc      + y.cv
                            cnst = x.cnst && y.cnst
                            for ip = 1:N
                                            cc_grad[ip] = scal * x.cv_grad[ip] + y.cc_grad[ip]
                                            cv_grad[ip] = scal * x.cc_grad[ip] + y.cv_grad[ip]
                            end
                            R[i] = MC{N}(cv,cc,IntervalType(lo,hi),SVector{N,Float64}(cv_grad), SVector{N,Float64}(cc_grad),cnst)
                    end
        end
        return R #Type Vector{MC}(undef, N)
end
