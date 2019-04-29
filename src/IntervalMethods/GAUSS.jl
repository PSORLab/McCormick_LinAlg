#This function estimates bounds on the determinate of an interval matrix
# by applying gaussian elimination
function gaussdet(A) #From A NEW CRITERION TO GUARANTEE THE FEASIBILITY OF THE INTERVAL GAUSSIAN ALGORITHM*, A. FROMMERf AND G. MAYER
      (m, n) = size(A);
      A_kplus1 = copy(A)
      A_k = copy(A)
      Intzero = Interval(0.0, 0.0)
      #Closer to paper's values
      for k = 1:(n-1)
            for i = 1:n #Consider Pivot Tightening?
                  for j = 1:n#Assume A is nxn
                        if 1<=i<=k && 1<=j<=n
                        elseif (k+1)<=i<=n && (k+1)<=j<=n
                              A_kplus1[i,j] = A_k[i,j] - (A_k[i,k]*A_k[k,j])/A_k[k,k]
                        else
                              A_kplus1[i,j] = Intzero
                        end
                  end
            end
            A_k = A_kplus1
      end
      det_k::Interval = A_k[1,1]
      for i = 2:n
            det_k *= A_k[i,i]
      end
      return det_k
end

function gaussddet(A) #Using dual extended IntArithmatic from T. Nirmala1, D. Datta2, H.S. Kushwaha3, K. Ganesan4 §
      (m, n) = size(A);
      Intzero = Interval(0.0, 0.0)
      Ad = copy(A)
      for i = 1:n #Assume A is nxn. No row swapping as of yet
            if Ad[i,i] != Intzero #As long as the current pivot isn't 0
                pivotval = Ad[i,i];
                for row=i:(m-1) #Row subtracts pivot from below to make it a true pivot column
                    if Ad[row+1,i] != 0#if the next col isnt 0
                         Ad[row+1,:] = intsub(Ad[row+1,:], intmult(intdiv(Ad[i,:], pivotval), Ad[row+1,i])); #Translates to newrow = oldrow - (pivot scaled to elimate pivot column)
                    end
                end
            end
      end
      detd::Interval = Ad[1,1]
      for i = 2:n
            detd = intmult(detd, Ad[i,i])
      end
      return detd
end

function gaussptdet(A)  #From A NEW CRITERION TO GUARANTEE THE FEASIBILITY OF THE INTERVAL GAUSSIAN ALGORITHM*, A. FROMMERf AND G. MAYER
      (m, n) = size(A); #Haven't found suitable test matrix for this yet
      A_kplus1 = copy(A)
      A_k = copy(A)
      Intzero = Interval(0.0, 0.0)
      #Closer to paper's values
      for k = 1:(n-1)
            if contains_zero(A[k,k])
                  A[k,k] = A[k,k] ∩ gaussptdet(A[1:k, 1:k])/gaussptdet(A[1:(k-1),1:(k-1)])
            end
            for i = 1:n #Consider Pivot Tightening?
                  for j = 1:n#Assume A is nxn
                        if 1<=i<=k && 1<=j<=n
                        elseif (k+1)<=i<=n && (k+1)<=j<=n
                              A_kplus1[i,j] = A_k[i,j] - (A_k[i,k]*A_k[k,j])/A_k[k,k]
                        else
                              A_kplus1[i,j] = Intzero
                        end
                  end
            end
            A_k = A_kplus1
      end
      det_k::Interval = A_k[1,1]
      for i = 2:n
            det_k *= A_k[i,i]
      end
      return det_k
end


#=
      function gaussdet(A)
      Intzero::Interval = Interval(0.0,0.0)
      (n,m) = size(A)
      Ad = copy(A_k)
      #Bounds are much tighter than expected
      for i = 1:n #Assume A is nxn. No row swapping as of yet
            if Ad[i,i] != Intzero #As long as the current pivot isn't 0
                pivotval = Ad[i,i];
                for row=i:(m-1) #Row subtracts pivot from below to make it a true pivot column
                    if Ad[row+1,i] != 0#if the next col isnt 0
                         Ad[row+1,:] = dualsub(Ad[row+1,:], Ad[row+1,i] * dualdiv(Ad[i,:], pivotval)); #Translates to newrow = oldrow - (pivot scaled to elimate pivot column)
                    end #Using dual from T. Nirmala1, D. Datta2, H.S. Kushwaha3, K. Ganesan4 §
                end
            end
      end
      detd::Interval = Ad[1,1]
      for i = 2:n
            detd *= Ad[i,i]
      end
       return detd
end
=#
