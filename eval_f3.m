function f= eval_f3(x,p,u,b, varargin)

num_bowties = p.NumBowties;

A = p.Area;
beta = p.Beta;
d = p.Distance;
phi = p.workFunction;

ROC = p.Radius;
taby = p.taby;
invC = p.invC;
CG = p.CG;
Ivec = zeros(2*num_bowties,1);

for i = 1:2*num_bowties
    if mod(i,2)==1
        Ivec(i) = - A*p.jnano*Jnano(phi,beta*(x(i)-x(i+1))/d,ROC,taby);
    else
        Ivec(i) = A*p.jnano*Jnano(phi,beta*(x(i-1)-x(i))/d,ROC,taby);
    end  
end

if length(varargin) > 1
    % this part is for linearization
    outputmode = varargin{2}; 
    if outputmode == 'onlyx0'
        f = CG*x + invC*Ivec;    
    else
        f = CG*x + invC*Ivec + invC*b*u;
    end
else
    f = CG*x + invC*Ivec + invC*b*u;
% CG*x
% invC*Ivec 
end
