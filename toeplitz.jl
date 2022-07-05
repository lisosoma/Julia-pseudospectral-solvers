function toeplitz(x::AbstractVector{T}) where T
   n = length(x)
   A = zeros(T, n, n)
   for i = 1:n
       for j = 1:n-i+1
           A[i,i+j-1] = x[j]
       end
       for j = n-i+2:n
           A[i, j-(n-i+1)] = x[j]
       end
   end
   return A
end