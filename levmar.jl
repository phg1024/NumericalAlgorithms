function jacobian(f, x)
  n = length(x);
  m = length(f(x));
  epsval = 1e-6;
  J = zeros(m, n);
  for i=1:n
    ivec = zeros(n, 1);
    ivec[i] = epsval;
    J[:, i]  = (f(x+ivec)-f(x))/epsval;
  end
  return J;
end

function levmar_dif(f, x0, kmax, opts)
  return levmar_der(f, jacobian, x0, kmax, opts);
end

function levmar_der(f, jac, x0, kmax, opts)
  k = 0;
  nv = 2;
  x = x0;
  n = length(x);
  m = length(f(x));

  eps1 = opts[1]; eps2 = opts[2]; eps3 = opts[3]; tau = opts[4];

  J = jac(f, x);
  A = J'*J; g = vec(J' * f(x));
  found = norm(g, Inf) <= eps1;
  mu = tau * maximum(diag(A));

  function F(x)
    return 0.5 * dot(f(x), f(x));
  end

  println("starting iterations...");
  while !found && k < kmax
    k += 1;
    h = vec((A + mu * eye(n, n)) \ (-g));

    if norm(h) <= eps2*(norm(x) + eps2)
      println("h too small");
      found = true
    else
      xnew = x + h;
      D = 0.5 * dot(h, mu*h-g);
      rho = (F(x) - F(xnew))/D;
      if rho > 0
        x = xnew;
        J = jac(f, x);
        A = J'*J; g = vec(J'*f(x));
        found = (norm(g, Inf) <= eps1) || (F(x) <= eps3);
        if found
          println("g too small");
        end
        mu = mu * max(1/3, 1 - (2*rho-1)^3);
        nu = 2;
      else
        mu = mu * nv; nv *= 2;
      end
    end
  end

  println("converged in ", k, " iterations. residue = ", F(x));
  return x;
end

function myequ(x)
  return [73.0 - x[1]^2, (4.0-x[2])^2];
end

function myequ2(x)
  return [exp(x[1]-2)+sin(x[1])*log(abs(x[1])), x[2]-2.0];
end

# Meyer's problem, minimum at (2.48, 6.18, 3.45)
function meyer(x)
  n = 16;
  f = zeros(n, 1);
  # The measurement vector
  mx = [34.78, 28.61, 23.65, 19.63,
        16.37, 13.72, 11.54, 9.744,
        8.261,  7.03, 6.005, 5.147,
        4.427,  3.82, 3.307, 2.872];
  for i=1:n
    ui = 0.45 + 0.05 * (i-1);
    f[i] = x[1]*exp(10.0*x[2]/(ui+x[3]) - 13.0) - mx[i];
  end
  return vec(f);
end

function jac_meyer(f, x)
  n = 16;
  J = zeros(n, 3);
  for i=1:n
    ui = 0.45 + 0.05 * (i-1);
    tmp = exp(10.0*x[2]/(ui+x[3]) - 13.0);
    J[i, 1] = tmp;
    J[i, 2] = 10.0*x[1]*tmp/(ui+x[3]);
    J[i, 3] = -10.0*x[1]*x[2]*tmp/((ui+x[3])*(ui+x[3]));
  end
  return J;
end

println(levmar_dif(myequ, [0.5, 0.5], 100, [1e-12, 1e-12, 1e-17, 1e-3]))
println(levmar_dif(myequ2, [0.5, 0.5], 100, [1e-8, 1e-8, 1e-17, 1e-3]))
println(levmar_dif(meyer, [8.85, 4.0, 2.5], 1000, [1e-17, 1e-17, 1e-17, 1e-3]))
println(levmar_der(meyer, jac_meyer, [8.85, 4.0, 2.5], 1000, [1e-17, 1e-17, 1e-10, 1e-3]))
