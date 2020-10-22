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
p.Cwire = 0.08;
p.Rwire = 1;
u.Wire1Bias = 10/p.Rwire;
u.Wire2Bias = 0;
u.jnano = 1;

v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;

x0 = zeros(4*p.NumBowties, 1);
% u0 = zeros(2*p.NumBowties, 1);
% u0(1:2:end) = v1/p.REmitter;
% u0(2:2:end) = v2/p.RCollector;

t_stop = 1;
t_start = 0;g
timestep = 0.0005;

% Implement forward euler
X = ForwardEulerNewest(x0,p,u,t_start,t_stop,timestep, @eval_f_new);

figure;
plot(X(2,:))
hold on;
plot(X(3,:))
plot(X(18, :))
plot(X(19, :))
plot(X(38, :))
plot(X(39, :))

legend('v left1','v right1', 'v left5','v right5', 'v left10','v right10')
title('voltages over time','FontSize', 8)
xlabel('time (ns)')
ylabel('node voltage (v)')
hold off;

