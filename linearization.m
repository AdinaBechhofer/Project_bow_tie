function [A,B] = linearization(f, x0, p, u)
epsilon = 0.001;
A = FiniteDifferenceJacobian(f, x0, p, u);
f_at_x0 = f(x0,p,u);
vE= u.vEmitter;
u.vEmitter = vE + epsilon;
Ju(:,1) = (f(x0,p,u) - f_at_x0)/epsilon;
u.vEmitter = vE;
vC = u.vCollector; 
u.vCollector = vC + epsilon;
Ju(:,2) = (f(x0,p,u) - f_at_x0)/epsilon;
u.vCollector = vC;
K0= f_at_x0 - A*x0 - Ju*[vE; vC];
%disp(K0)
B = [K0,Ju];
%disp(A*x0- f(x0,p,u))

end
