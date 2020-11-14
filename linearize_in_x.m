function [A,B] = linearize_in_x(f, x0, t, p, u, b)
%epsilon = 0.001;
A = FiniteDifferenceJacobian_t(f, x0, t, p,u,b);
f_at_x0 = f(x0,t,p,u,b);
%vE= u.vEmitter;
%u.vEmitter = vE + epsilon;
%Ju(:,1) = (f(x0,p,u) - f_at_x0)/epsilon;
%u.vEmitter = vE;
%vC = u.vCollector; 
%u.vCollector = vC + epsilon;
%Ju(:,2) = (f(x0,p,u) - f_at_x0)/epsilon;
%u.vCollector = vC;
K0= f_at_x0 - A*x0;
%disp(K0)
B = [K0,b];
%disp(A*x0- f(x0,p,u))

end
