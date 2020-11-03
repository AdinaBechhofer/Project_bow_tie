function out = f_simple(x, p, u)
out = x.*sin(x) +2*u.vEmitter^2 -0.2*sqrt(u.vCollector);
end 