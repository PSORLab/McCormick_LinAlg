#This function estimates bounds on the determinate of an interval matrix
# by applying gaussian elimination
using IntervalArithmetic
function dual(Inter::Interval) # E. Kaucher, Interval analysis in the extended interval space IR, Computing, Suppl., 2 (1980), 33-49, doi: 10.1007/978-3-7091-8577-3
      return Interval(Inter.hi, Inter.lo)
end
function dualsub(Inter1, Inter2)
      return Interval(Inter1.lo - Inter2.lo, Inter1.hi-Inter2.hi)
end
function dualdiv(Inter1, Inter2)
      return Inter1 / dual(Inter2)
end
function test()
Intzero::Interval = Interval(0.0,0.0)
Am = [1. 2. 3.
      4. 6. 7.
      5. 9. 8.]
r1 = .1
r2 = .01
A_ = map(x->IntervalType(x-r1, x+r1), Am) #Interval Matrix
HULL = [[4.060, 14.880] [8.465, 9.545]]
#Begin upper-trianglulating the matrix
(m, n) = size(A_);
detf::Float64 = 1 #Accumulated determinatn factor
A_kplus1 = copy(A_)
A_k = copy(A_)
for k = 1:(n-1) #From A NEW CRITERION TO GUARANTEE THE FEASIBILITY OF THE INTERVAL GAUSSIAN ALGORITHM*, A. FROMMERf AND G. MAYER
      for i = 1:n
            for j = 1:n#Assume A is nxn
                  if 1<=i<=k && 1<=j<=n
                  elseif (k+1)<=i<=(k+n) && (k+1)<=j<=(k+n)
                        A_kplus1[i,j] -= (A_[i,k]*A_[k,j])/A_[k,k]
                  else
                        A_kplus1[i,j] = 0
                  end
            end
      end
      A_ = A_kplus1
end

Ad = copy(A_)
for i = 1:n #Assume A is nxn. No row swapping as of yet
      if Ad[i,i] != Intzero #As long as the current pivot isn't 0
          pivotval = Ad[i,i];
          for row=i:(m-1) #Row subtracts pivot from below to make it a true pivot column
              if Ad[row+1,i] != 0#if the next col isnt 0
                   Ad[row+1,:] = dualsub(Ad[row+1,:], Ad[row+1,i] * dualdiv(Ad[i,:] / pivotval)); #Translates to newrow = oldrow - (pivot scaled to elimate pivot column)
              end #Using dual from T. Nirmala1, D. Datta2, H.S. Kushwaha3, K. Ganesan4 §
          end
      end
end
#=
for i=1:m #checks every row
    LV = abs(A[i, i]); #largest value set to first item as default
    row = i; #Default for row index to be swapped following the for loop
    #PARTIAL PIVOTING
    for j=(i+1):m#Finds largest valued pivot in current colomn
        if !(LV >= abs(A[j, i]))
            LV = abs(A[j, i]);
            row = j;
        end
    end
    if row != i #Applys the row swap so long as there was a row with larger pivot value than current
        placeholder = A[i,:];
        A[i,:] = A[row,:];
        A[row,:] = placeholder;
    end
    #ELIMINATION
    if A[i,i] != 0 #As long as the current pivot isn't 0, which would mean the entire column is zero's
        pivotval = A[i,i];
        for row=i:(m-1) #Row subtracts pivot from below to make it a true pivot column
            if A[row+1,i] != 0#if the next col isnt 0
                 A[row+1,:] = A[row+1,:] - A[i,:] * (A[row+1,i] / pivotval); #Translates to newrow = oldrow - (pivot scaled to elimate pivot column)
            end
        end
        #A checkpoint after elim
    end
end
=#
det_::Interval = A_[1,1]
detd::Interval = Ad[1,1]
for i = 2:n
      det_ *= A_[i,i]
      detd *= Ad[i,i]
end
 return A_, det, Ad, detd
end
A_, det, Ad, detd = test()
