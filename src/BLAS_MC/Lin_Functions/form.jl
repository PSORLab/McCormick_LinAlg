#Functions to take matrices from full storage to X - style storage

function band(AF::Array{MC{N},2},m,n,ku,kl) where N #undef makes messy matrix with e-300 instead of 0's
    AB::Array{MC{N},2} = Array{MC{N}}(undef, m,n) #looks messy with undef values isntead of 0's
    for j = 1:n
        k = ku + 1 - j
        for i = max(1,j-ku):min(m,j+kl)
            AB[k+i,j] = AF[i,j]
        end
    end
    return AB
end
function sbandu(AF::Array{MC{N},2},m,n,k) where N
    AB::Array{MC{N},2} = Array{MC{N}}(undef, m,n) #looks messy with undef values isntead of 0's
    for j = 1:n
        M = k + 1 - j
        for i = max(1,j-k):j
            AB[M+i,j] = AF[i,j]
        end
    end
    return AB
end
function sbandl(AF::Array{MC{N},2},m,n,k) where N
    AB::Array{MC{N},2} = Array{MC{N}}(undef, m,n) #looks messy with undef values isntead of 0's
    for j = 1:n
        M = 1 - j
        for i = j:min(n,j+k)
            AB[M+i,j] = AF[i,j]
        end
    end
    return AB
end
