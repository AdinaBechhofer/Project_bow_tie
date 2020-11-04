clear all;
close all;

p.NumBowties = 100;
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
% number of rows
p.row = 10;
% number of columns 
p.col = 10;
p.Ccoupling = 0.03;
%u.jnano = 1;

v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;

C = zeros(2*p.NumBowties);
G = zeros(2*p.NumBowties);
% Stamp C and invert; stamp G
for i = 1:2*p.NumBowties
    C(i,i) = p.CemitterCollector+p.Cparasitic;
    if mod(i,2) == 1
        C(i,i+1)= C(i,i+1)-p.CemitterCollector;
        G(i,i) = G(i,i) +1./(p.REmitter);
        if  i> 2*p.row
            C(i,i) = C(i,i) + p.Ccoupling;
            C(i,i-2*p.row+1) = C(i,i-2*p.row+1) -p.Ccoupling;
        end 
    else
        C(i,i-1)= C(i,i-1) -p.CemitterCollector;
        G(i,i) = 1./(p.RCollector);
        if i<= 2*p.row*(p.col-1)
            C(i,i) = C(i,i) + p.Ccoupling;
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) - p.Ccoupling;
        end
    end
end
figure;
spy(C)
issymmetric(C)
mean(sum(C,2))
issymmetric(G)
p.invC = inv(C);
p.CG = p.invC*G;

x0 = zeros(2*p.NumBowties, 1);

t_stop = 0.1;
t_start = 0;
timestep = 0.000005;

% Implement forward euler
X = ForwardEulerNew(x0,p,u,t_start,t_stop,timestep);

figure;git
plot(X(5,:))
hold on;
plot(X(6,:))
plot(X(5+4*p.col,:))
plot(X(6+4*p.col,:))
plot(X(5+2*(p.row-1)*p.col, :))
plot(X(6+2*(p.row-1)*p.col, :))
legend('vodd1','veven1','vodd2','veven2','vodd3','veven3')
title('voltages over time','FontSize', 8)
hold off;

