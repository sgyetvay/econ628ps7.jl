module econ628ps7

using ForwaDiff, LinearAlgebra

function newton(f, f_prime, x_0; tolerance = 1E-7, maxiter = 1000)
    #Setup the algorithm
    x_old = x_0
    error = Inf
    iter = 1
    while error > tolerance && iter <= maxiter
        x_new = x_old - (f(x_old)/f_prime(x_old))
        error = sqrt((x_new - x_old)^2)
        x_old = x_new
        iter = iter + 1
    end
    if iter == maxiter + 1
        return nothing
    else
        return x_old
end
end

f(x) = (x-1)^3
f_prime(x) = 3*(x-1)^2

function newton(f, x_0; tolerance = 1E-7, maxiter = 1000)
f_prime = x -> ForwardDiff.derivative(f, x)
return newton(f,f_prime,x_0;tolerance = tolerance, maxiter = maxiter)
end

export newton, f, f_prime
end
