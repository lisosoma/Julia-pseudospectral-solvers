function DN_generator(N)
    h = 2 * pi / N
    M = [mod(j, N) != 0 ? (-1) ^ (j+1) * cot(j * h / 2) / 2 : 0 for j = 0:(N - 1)]
    DN = transpose.(toeplitz(M))
    return DN
end