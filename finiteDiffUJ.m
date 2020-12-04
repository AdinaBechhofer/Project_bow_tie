function [A,B] = finiteDiffUJ(f, x0, p, u)
epsilon = 0.01;
A = FiniteDifferenceJacobian(@eval_f2, x0, p, u);
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
B = [K0,Ju];
end