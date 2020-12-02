function [F, JF] = FJFTrap(x, p, u, b, varargin)
t = varargin{1};
gamma = varargin{2};
dt = varargin{3};
fJfhand = varargin{4};

[f, Jf] = fJfhand(x,p, u, b, t);
F = x - (dt/2)*f - gamma;
JF = eye(length(Jf)) - (dt/2)*Jf;
%disp(cond(Jf))
end