p.NumBowties = 50;
p.REmitter = 1; % 1 ohm 
p.RCollector = 1; % 1 ohm 
p.Area = 0.0001; % need to actually get a reasonable estimate  
p.Beta = 0.9; % enhancement factor 
p.Distance = 10e-7; % 10 nm
p.workFunction = 8.49e-19; % work function of gold 
p.CemitterCollector = 1e-9; % nano farad
p.Cparasitic = 1e-10; % 0.1 nano farad
p.Radius = 1e-7; % 1 nm?


x0 = zeros(2*p.NumBowties, 1);
u0 = zeros(2*p.NumBowties, 1);
u0(1:2:end) = 5;
u0(2:2:end) = 0.5;

t_stop = 400;
t_start = 1;
timestep = 1;

% Implement forward euler
X = ForwardEulerNew(x0,p,u0,t_start,t_stop,timestep);

figure;
plot(X(5,:))
hold on;
plot(X(6,:))

legend('vodd','veven')
title('voltages over time','FontSize', 8)
hold off;

