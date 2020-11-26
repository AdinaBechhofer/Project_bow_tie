function [A,B] = linearization(f,x0,t,p,u,b)
epsilon = 0.001;
A = FiniteDifferenceJacobian_t(f, x0,t, p, u,b);
f_at_x0 = f(x0,t,p,u,b);
Ju = zeros(length(x0), length(u));
for k=1:length(u)
    ek = zeros(size(u));
    ek(k) = 1;
    Ju(:,k) = (f(x0,t,p,u+epsilon*ek,b) - f_at_x0)/epsilon;
end
K0= f_at_x0 - A*x0 - Ju*u;
B = [K0,Ju];
%disp(A*x0- f(x0,p,u))

end
