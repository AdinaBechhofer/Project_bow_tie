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

x0 = zeros(2*p.NumBowties, 1);

t_stop = 1;
t_start = 0;
timestep = 0.05;

% Implement forward euler
X = ForwardEulerNew(x0,p,u,t_start,t_stop,timestep);

x0 = X(:,1);
x0 = X(:,5);
x0 = X(:,20);

FDJ = FiniteDifferenceJacobian(@eval_f, x0, p, u, 0.01);
disp(cond(FDJ))