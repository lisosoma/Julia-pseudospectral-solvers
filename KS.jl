import Pkg
Pkg.add("Plots")
Pkg.add("FFTW")
using Plots
using FFTW

#----------------------------------------------------#
# functions for numerical method solving nonlinear pde

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


function DN_generator(N, l)
    h = l * pi / N
    M = [mod(j, N) != 0 ? (-1) ^ (j+1) * 
         cot(j * h / l) / l : 0 for j = 0:(N - 1)]
    DN = transpose.(toeplitz(M))
    return DN
end

#----------------------------------------------------#
# function for right part of KS equation

function phi(F, l)
    return [F[k] * (2 * pi * k / (l * pi))^2 -  F[k] * (2 * pi * k / (l * pi))^4 for k=1:length(F)]
end


#----------------------------------------------------#
# solution

N = 200 # number of grid points
l = 30
DN = DN_generator(N, l)

x = zeros(N) # points of grid
x = [l * pi * j / N for j = 1:N]

T = 200 # time
h = ((l * pi / N) ^ 2) / 300

u0 = [[(cos(l*pi*x[i] / 2)) for i=1:N]]
for j = 1:round(Int, T / h - 1)
    u0_hat = rfft(last(u0))
    u1_hat = u0_hat + h * phi(u0_hat, l)
    u1 = -h * (last(u0) .* (DN * last(u0))) + irfft(u1_hat, N)
    push!(u0, u1)
end

M4 = reduce(vcat,transpose.(u0))
heatmap(x, t, M4, dpi=300)
#savefig("heatmap_KS7.png")
