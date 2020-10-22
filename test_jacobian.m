clear all;
close all;

p.NumBowties = 100;
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
p.Col = 10;
p.Row = 10;
p.Cwire = 0.08*1e5;
p.Rwire = 1*1e-4;
u.Wire1Bias = 10/p.Rwire;
u.Wire2Bias = 0;
u.jnano = 1;

x0 = zeros(4*p.NumBowties, 1);

t_stop = 5;
t_start = 0;
timestep = 0.0005;

% Implement forward euler
X = ForwardEulerNewest(x0,p,u,t_start,t_stop,timestep, @eval_f_new);

for i = [1,1000,5000]
    
    x0 = X(:,i);
    %x0 = X(:,1000);
    %x0 = X(:,5000);

    FDJ = FiniteDifferenceJacobian(@eval_f_new, x0, p, u, 0.01);
    spy(FDJ)
    disp(cond(FDJ))
end