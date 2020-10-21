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
% u0 = zeros(2*p.NumBowties, 1);
% u0(1:2:end) = v1/p.REmitter;
% u0(2:2:end) = v2/p.RCollector;

t_stop = 1;
t_start = 0;
timestep = 0.05;

% Implement forward euler
X = ForwardEulerNew(x0,p,u,t_start,t_stop,timestep);

figure;
plot(X(5,:))
hold on;
plot(X(6,:))

legend('vodd','veven')
title('voltages over time','FontSize', 8)
hold off;

