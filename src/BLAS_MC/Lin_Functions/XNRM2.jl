#This code is based on the BLAS implementation of SDOT
#At this point is a crude reimplementation of xDOT(X,Y) for X = Y

#Haven't included strided vectors yet, more important on matrix functions
#::SVector{n, MC{N}}
function XNRM2(X::SVector) #x,y E(Vector(MC{N})) where N <: Integer
    n::Int = length(X)
    N::Int = length(X[1].cv_grad)
    cum_cc::Float64 = 0.0
    cum_cv::Float64 = 0.0
    cum_hi::Float64 = 0.0
    cum_lo::Float64 = 0.0
    cum_ccgrad::Vector{Float64} = Vector{Float64}(undef,N)
    cum_ccgrad::Vector{Float64} .= 0
    cum_cvgrad::Vector{Float64} = Vector{Float64}(undef,N)
    cum_cvgrad::Vector{Float64} .= 0
    cum_cnst::Bool = 1

    m::Integer = mod(n,5)
    if m != 0
        for i in 1:m
            #global cum_cv, cum_cc, cum_cv, cum_hi, cum_lo, cum_cvgrad, cum_ccgrad, cum_const
        #= #If multtemp returns list of values
           temp::SVector{6,Any} = multtemp(X[i], Y[i])
            cum_cv += temp[1]
            cum_cc += temp[2]
            cum_hi += temp[3].hi
            cum_lo += temp[3].lo
            cum_cvgrad += temp[4]
            cum_ccgrad += temp[5]
            cum_const = (cum_const && temp[6])
            =#
            #Still passing MC{N}'s
            temp::MC = X[i]^2 #Mult function needs more integration here
            #temp = multtemp(X[i], Y[i])
            cum_cv += temp.cv
            cum_cc += temp.cc
            cum_hi += temp.Intv.hi
            cum_lo += temp.Intv.lo
            cum_cvgrad += temp.cv_grad #Vector += SVector
            cum_ccgrad += temp.cc_grad
            cum_cnst = (cum_cnst && temp.cnst)
        end
    end
        #now continue in series of 5 up to i = n-m+1 , m from mod(n,m)
        for i in (m+1):5:n
            #global cum_cv, cum_cc, cum_cv, cum_hi, cum_lo, cum_cvgrad, cum_ccgrad, cum_const
        #= #If multtemp returns list of values
           temp1::SVector{6,Any} = multtemp(X[i], Y[i])
           temp2::SVector{6,Any} = multtemp(X[i+1], Y[i+1])
           temp3::SVector{6,Any} = multtemp(X[i+2], Y[i+2])
           temp4::SVector{6,Any} = multtemp(X[i+3], Y[i+3])
           temp5::SVector{6,Any} = multtemp(X[i+4], Y[i+4])
           cum_cv += temp1[1] +temp2[1] +temp3[1] +temp4[1] +temp5[1]
           cum_cc += temp1[2] +temp2[2] +temp3[2] +temp4[2] +temp5[2]
           cum_hi += temp1[3].hi +temp2[3].hi +temp3[3].hi +temp4[3].hi +temp5[3].hi
           cum_lo += temp1[3].lo +temp2[3].lo +temp3[3].lo +temp4[3].lo +temp5[3].lo
           cum_cvgrad += temp1[4] +temp2[4] +temp3[4] +temp4[4] +temp5[4]
           cum_ccgrad += temp1[5] +temp2[5] +temp3[5] +temp4[5] +temp5[5]
           cum_const = (cum_const && temp1[6] && temp2[6] && temp3[6] && temp4[6] && temp5[6])
           =#
           #Still passing MC{N}'s
           temp1::MC = X[i]^2
           temp2::MC = X[i+1]^2
           temp3::MC = X[i+2]^2
           temp4::MC = X[i+3]^2
           temp5::MC = X[i+4]^2
           cum_cv += temp1.cv +temp2.cv +temp3.cv +temp4.cv +temp5.cv
           cum_cc += temp1.cc +temp2.cc +temp3.cc +temp4.cc +temp5.cc
           cum_hi += temp1.Intv.hi +temp2.Intv.hi +temp3.Intv.hi +temp4.Intv.hi +temp5.Intv.hi
           cum_lo += temp1.Intv.lo +temp2.Intv.lo +temp3.Intv.lo +temp4.Intv.lo +temp5.Intv.lo
           cum_cvgrad += temp1.cv_grad +temp2.cv_grad +temp3.cv_grad +temp4.cv_grad +temp5.cv_grad
           cum_ccgrad += temp1.cc_grad +temp2.cc_grad +temp3.cc_grad +temp4.cc_grad +temp5.cc_grad
           cum_cnst = (cum_cnst && temp1.cnst && temp2.cnst && temp3.cnst && temp4.cnst && temp5.cnst)
       end

    result::MC = MC{N}(cum_cv, cum_cc, IntervalType(cum_lo, cum_hi), SVector{N,Float64}(cum_cvgrad), SVector{N,Float64}(cum_ccgrad), cum_cnst)#MCCormick Object
    result = sqrt(result)  #May not be opt way to call this
    return result
end
