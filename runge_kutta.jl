function runge_kutta(w, tn, f, h)
    
    for j = 1:round(Int, tn / h - 1)
        k1 = h * f(w[j])
        k2 = h * f((w[j] + k1) / 2)
        k3 = h * f((w[j] + k2) / 2)
        k4 = h * f(w[j] + k3)
        
        push!(w, w[j] + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
    end
    
    return w
end