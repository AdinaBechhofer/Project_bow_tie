function [A,B] = linearization(f,x0,p,u,b, varargin)
epsilon = 0.001;
if ~isempty(varargin)
    t = varargin{1}
    A = FiniteDifferenceJacobian_t(f, x0,p, u,b, t);
    f_at_x0 = f(x0,p,u,b,t);
Ju = zeros(length(x0), length(u));
for k=1:length(u)
    ek = zeros(size(u));
    ek(k) = 1;
    Ju(:,k) = (f(x0,p,u+epsilon*ek,b,t) - f_at_x0)/epsilon;
end
else 
    A = FiniteDifferenceJacobian_t(f, x0,p, u,b);
    f_at_x0 = f(x0,p,u,b);
Ju = zeros(length(x0), length(u));
for k=1:length(u)
    ek = zeros(size(u));
    ek(k) = 1;
    Ju(:,k) = (f(x0,p,u+epsilon*ek,b) - f_at_x0)/epsilon;
end 
end 
K0= f_at_x0 - A*x0 - Ju*u;
B = [K0,Ju];
%disp(A*x0- f(x0,p,u))

end
