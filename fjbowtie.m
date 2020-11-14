function [f, J] = fjbowtie(x,t,p,u,b)
T = u.period; %timestep approx 0.001
u_p = u.amp*cos(2*pi/T *t);

uEo = u.vEmitter;
uCo = u.vCollector;
u.vEmitter = uEo * (0 + u_p);
u.vCollector = uCo * (0+ u_p);
f = eval_f3(x,t,p,u,b);
J = FiniteDifferenceJacobian_t(@eval_f3, x, t, p, u, b);
u.vEmitter = uEo;
u.vCollector = uCo;
