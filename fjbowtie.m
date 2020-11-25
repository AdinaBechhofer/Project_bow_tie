function [f, J] = fjbowtie(x,t,p,u,b)
T = u.period; %timestep approx 0.001

u.sinuE = u.ampE*cos(2*pi/T *t + u.phaseE);
u.sinuC = u.ampC*cos(2*pi/T *t + u.phaseC);

%uEo = u.vEmitter;
%uCo = u.vCollector;
%u.vEmitter = uEo * (1 + u_p);
%u.vCollector = uCo * (1 + u_p);
f = eval_f3(x,t,p,u,b);
J = FiniteDifferenceJacobian_t(@eval_f3, x, t, p, u, b);
%u.vEmitter = uEo;
%u.vCollector = uCo;
