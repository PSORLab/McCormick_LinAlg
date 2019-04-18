#This function estimates bounds on the determinate of an interval matrix
# by applying gaussian elimination
function gaussdet(A)
(m, n) = size(A);
#Best Method
A_kplus1 = copy(A)
A_k = copy(A)
#Closer to paper's values
for k = 1:(n-1) #From A NEW CRITERION TO GUARANTEE THE FEASIBILITY OF THE INTERVAL GAUSSIAN ALGORITHM*, A. FROMMERf AND G. MAYER
      for i = 1:n
            for j = 1:n#Assume A is nxn
                  if 1<=i<=k && 1<=j<=n
                  elseif (k+1)<=i<=(k+n) && (k+1)<=j<=(k+n)
                        A_kplus1[i,j] -= (A_k[i,k]*A_k[k,j])/A_k[k,k]
                  else
                        A_kplus1[i,j] = 0
                  end
            end
      end
      A_k = A_kplus1
end

#Bounds are much tighter than expected
Ad::Array{Interval{Float64},2} = copy(A)
pivotval::Interval{Float64} = Interval(0.0, 0.0)
for i = 1:n #Assume A is nxn. No row swapping as of yet
      if Ad[i,i] != Intzero #As long as the current pivot isn't 0
          pivotval = Ad[i,i];
          for row=i:(m-1) #Row subtracts pivot from below to make it a true pivot column
              if Ad[row+1,i] != 0#if the next col isnt 0
                    #println(i, " ", row)
                   Ad[row+1,:] = dualsub(Ad[row+1,:], Ad[row+1,i] * dualdiv(Ad[i,:], pivotval)); #Translates to newrow = oldrow - (pivot scaled to elimate pivot column)
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
#println()
#println()
#display(Ad)
#display(A_k)
detd::Interval = Ad[1,1]
det_k::Interval = A_k[1,1]
for i = 2:n
      detd *= Ad[i,i]
      det_k *= A_k[i,i]
end
 return detd, det_k
end
