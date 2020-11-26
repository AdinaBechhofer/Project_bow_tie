function out = f_simple(x, t,p,u,b)
out = x.*sin(x) +2*u(1)^2 -0.2*sqrt(u(2));
end 