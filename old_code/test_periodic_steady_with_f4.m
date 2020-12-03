clear all;
close all;
tic;
p.NumBowties = 10;
% number of rows
p.row = 5;
% number of columns 
p.col = p.NumBowties/p.row;

p.REmitter = 1; % 1 ohm 
p.RCollector = 1; % 1 ohm 
%p.Area = 0.0001; % need to actually get a reasonable estimate  
p.Area = 5e-18; % 5 nm^2
%p.Beta = 0.9; % enhancement factor 
p.Beta = 25; % enhancement factor 
p.Distance = 10; % 10 nm
p.workFunction = 5.1; % work function of gold 
p.CemitterCollector = 2; % nano farad
p.Cparasitic = 0.05; % 0.1 nano farad
% p.Radius = 1; % 1 nm?
p.Radius = 10; % 1 nm
p.taby = csvread('rspa20140811supp3.csv');    

p.Ccoupling = 0.03;
%u.jnano = 1;


v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;
unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);

x0 = ones(2*p.NumBowties, 1);

t_stop = 0.5;
t_start = 0;
timestep = 0.001;

tvec = t_start:timestep:t_stop;
u.period = 20*timestep;

% %u.amp = 0;
% %X0 = ForwardEuler_t(@fjbowtie,x0,p,u,tvec,b);
% %X0 = TrapMethod(x0,@(xp,tp)fjbowtie(xp,tp,p,u,b), tvec);
u.amp = 1;
% %X = ForwardEuler_t(@fjbowtie,x0,p,u,tvec,b);
X = TrapMethod(x0,@(xp,tp)fjbowtie(xp,tp,p,u,b), tvec);


figure(1);
plot(tvec,X(1,:),'--')
hold on;
plot(tvec,X(2,:))

x0_new = newtonS(@(xp,tp)fjbowtie(xp,tp,p,u,b),x0,u,1);
xT = TrapMethod(x0_new,@(xp,tp)fjbowtie(xp,tp,p,u,b),tvec);

plot(tvec, xT(1,:),'--',tvec, xT(2,:))
legend('vodd-before','veven-before','vodd-periodic','veven-periodic')

hold off;
% plot(tvec,X(6,:), tvec,xT(6,:))
% hold off
% ylabel('v(t)')
% xlabel('t')
% legend('vodd-before','vodd-steady-state','veven-before','veven-steady-state')


% plot(X(5+4*p.col,:))
% plot(X(6+4*p.col,:))
% plot(X(5+2*(p.row-1)*p.col, :))
% plot(X(6+2*(p.row-1)*p.col, :))
% legend('vodd1','veven1','vodd2','veven2','vodd3','veven3')
% title('voltages over time','FontSize', 8)
% hold off;




end_time = toc