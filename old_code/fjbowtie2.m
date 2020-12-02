function [f, J] = fjbowtie2(x,t,p,u,b)
f = eval_f4(x,t,p,u,b); 
J = FiniteDifferenceJacobian_t(@eval_f4, x, t, p, u, b);

