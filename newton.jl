function derivative(f, x)
  eps = 1e-12
  return (f(x+eps)-f(x-eps))/(2.0*eps)
end

function myequ(x)
  return 73.0 - x^2
end

function myequ2(x)
  return exp(x-2)+sin(x)*log(abs(x))
end

function newtonsolve(equ, thres, maxiter)
  x = rand()
  iters = 0
  while abs(equ(x)) > thres && iters < maxiter
    println(x, '\t', equ(x), '\t', derivative(equ, x))
    x = x - equ(x)/derivative(equ, x)
    iters += 1
  end
  println("error = ", abs(equ(x)))
  println("converged in ", iters, " iters")
  return x
end

println(newtonsolve(myequ, 1e-12, 100))
println(newtonsolve(myequ2, 1e-12, 100))
