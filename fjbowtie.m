function [f, J] = fjbowtie(x,p,u,b, varargin)
%T = p.period; %timestep approx 0.001

%u.sinuE = p.ampE*cos(2*pi/T *t + p.phaseE);
%u.sinuC = p.ampC*cos(2*pi/T *t + p.phaseC);

%uEo = u.vEmitter;
%uCo = u.vCollector;
%u.vEmitter = uEo * (1 + u_p);
%u.vCollector = uCo * (1 + u_p);
if ~isempty(varargin)
    t = varargin{1};
    f = eval_f3(x,p,u,b, t);
    J = FiniteDifferenceJacobian_t(@eval_f3, x, p, u, b, t);
else 
    f = eval_f3(x,p,u,b);
    J = FiniteDifferenceJacobian_t(@eval_f3, x, p, u, b);
end 
%u.vEmitter = uEo;
%u.vCollector = uCo;
