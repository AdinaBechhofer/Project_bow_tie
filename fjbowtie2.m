function [f, J] = fjbowtie2(x,t,p,u,b)
%Su_p = u.amp*cos(2*pi/T *t);

% uEo = u.vEmitter;
% uCo = u.vCollector;
% u.vEmitter = uEo * (1 + u_p);
% u.vCollector = uCo * (0+ u_p);
f = eval_f4(x,t,p,u,b);
J = FiniteDifferenceJacobian_t(@eval_f4, x, t, p, u, b);

