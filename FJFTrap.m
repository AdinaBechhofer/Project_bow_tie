function [F JF] = FJFTrap(fJfhand, x, t, gamma,dt)

[f, Jf] = fJfhand(x,t);
F = x - (dt/2)*f - gamma;
JF = eye(length(Jf)) - (dt/2)*Jf;
%disp(cond(Jf))
end