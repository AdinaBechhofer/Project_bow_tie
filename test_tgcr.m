clear all;
close all;

p.NumBowties = 50;
p.REmitter = 1; % 1 ohm 
p.RCollector = 1; % 1 ohm 
p.Area = 0.0001; % need to actually get a reasonable estimate  
p.Beta = 0.9; % enhancement factor 
p.Distance = 10; % 10 nm
p.workFunction = 5.1; % work function of gold 
p.CemitterCollector = 0.1; % nano farad
p.Cparasitic = 0.08; % 0.1 nano farad
p.Radius = 1; % 1 nm?
p.taby = csvread('rspa20140811supp3.csv');    

u.jnano = 1;

v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;

% Stamp C and invert; stamp G
C = zeros(2*p.NumBowties);
G = zeros(2*p.NumBowties); % if G is non-linear need to stamp inside eval_f 
for i = 1:2*p.NumBowties
    C(i,i) = p.CemitterCollector+p.Cparasitic;
    if mod(i,2) == 1
        C(i,i+1)= -p.CemitterCollector;
        G(i,i) = -1./(p.REmitter);
    else
        C(i,i-1)= -p.CemitterCollector;
        G(i,i) = -1./(p.RCollector);
    end
end
invC = inv(C);
p.CG = invC*G;

x0 = zeros(2*p.NumBowties, 1);

t_stop = 5;
t_start = 0;
timestep = 0.01;

% Implement forward euler
X = ForwardEulerNew(x0,p,u,t_start,t_stop,timestep);

x0 = X(:,100);
%x0 = X(:,5);
%x0 = X(:,500);

f= -eval_f2(x0, p, u);
J = FiniteDifferenceJacobian(@eval_f2, x0, p, u, 0.01);
delX = J\f;

tol = 1e-5;
[x,r] = tgcr(J,f,tol,100);

maxeval = max(abs(eig(J)))
mineval = min(abs(eig(J)))

% Chebyshev bounds: 
A = (sqrt(maxeval/mineval)-1)/(sqrt(maxeval/mineval)+1);
log(tol*2)/log(A)