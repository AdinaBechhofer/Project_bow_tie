function [F, JF] = FJFTrap(x, p, u, b, varargin)
%disp(varargin)
input_args = varargin{1};
t = input_args{1};
gamma = input_args{2};
dt = input_args{3};
fJfhand = input_args{4};

[f, Jf] = fJfhand(x,p, u, b, t);
F = x - (dt/2)*f - gamma;
JF = eye(length(Jf)) - (dt/2)*Jf;
%disp(cond(Jf))
end