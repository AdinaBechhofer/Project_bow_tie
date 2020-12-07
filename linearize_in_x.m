function [A,B] = linearize_in_x(f, x0, p, u, b, varargin)
if ~isempty(varargin)
    t = varargin{1};
    A = FiniteDifferenceJacobian_t(f, x0,p,u,b,t);
    f_at_x0 = f(x0,p,u,b,t,'nonlinearpart');
else
% A = FiniteDifferenceJacobian_t(f, x0,p,u,b);
% f_at_x0 = f(x0,p,u,b);
    error("add nonlinearpart flag")
end

K0= f_at_x0 - A*x0;
B = [K0,p.invC*b]; % B is invC*b in the formulation

end
