clear all;
close all;

p.NumBowties = 100;
% number of rows
p.row = 10;
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
p.jnano = 1;


% Stamp C and invert; stamp G
C = zeros(2*p.NumBowties);
G = zeros(2*p.NumBowties); % this is for the linear part. Non-linear part is stamped in eval_f 
for i = 1:2*p.NumBowties
    C(i,i) = p.CemitterCollector+p.Cparasitic;
    if mod(i,2) == 1
        C(i,i+1)= C(i,i+1)-p.CemitterCollector;
        G(i,i) = G(i,i) -1./(p.REmitter);
        if  i> 2*p.row
            C(i,i) = C(i,i) + p.Ccoupling;
            C(i,i-2*p.row+1) = C(i,i-2*p.row+1) -p.Ccoupling;
        end 
    else
        C(i,i-1)= -p.CemitterCollector;
        G(i,i) = -1./(p.RCollector);
        if i<= 2*p.row*(p.col-1)
            G(i,i) = G(i,i) + p.Ccoupling;
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) - p.Ccoupling;
        end
    end
end
% figure;
% spy(C)
p.invC = inv(C);
p.CG = p.invC*G;

v1 = 5;
v2 = 0;
period =1.5;
amplitude = 0.5;
u = [v1/p.REmitter; v2/p.RCollector];
unitb = [1, 0, 1, 0; 0, 1, 0, 1];
b = repmat(unitb,p.NumBowties,1);
b(50:end, 3) = 0;
b(1:end-50, 4) = 0;

x0 = zeros(2*p.NumBowties, 1);
tic;
t_stop = 20;
t_start = 0;
timestep = 0.0005;
tvec_normal = t_start:timestep:t_stop;
U = [repmat(u,1,length(tvec_normal)); amplitude*cos(2*pi*tvec_normal/period); amplitude/2*cos(pi*tvec_normal/period +0.03)];
%U = [repmat(u,1,length(tvec_normal)); zeros(2, length(tvec_normal))];
X_normal = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec_normal);
%X = TrapMethod(x0,p,U,b,@fjbowtie, tvec);
time_euler = toc

tic;
tvec = 0:0.5:period ;
x0 = repmat([4.8; 0.2], p.NumBowties, 1);
U = [repmat(u,1,length(tvec)); amplitude*cos(2*pi*tvec/period); amplitude/2*cos(2*pi*tvec/period +0.03)];
x0_new = newtonS(x0,p,U,b,@fjbowtie,tvec);
time_shooting = toc 

tic;
tvec = t_start:timestep:t_stop;
U = [repmat(u,1,length(tvec)); amplitude*cos(2*pi*tvec/period); amplitude/2*cos(2*pi*tvec/period +0.03)];
X_ss = ForwardEuler_t(@eval_f3,x0_new,p,U,b,tvec);
time_post_euler = toc

figure;
subplot(2, 1, 1)
plot(tvec_normal, X_normal(1, :), tvec_normal, X_normal(2, :), 'linewidth', 1.2)
hold on
plot(tvec, X_ss(1,:),'--',tvec, X_ss(2,:), '--', 'linewidth', 1.2)
legend('v1-normal','v2-normal','v1-ss','v2-ss')
hold off;
ylabel('Voltage (V)')
xlabel('time (ns)')
subplot(2, 1, 2)
plot(tvec_normal, X_normal(end-1, :), tvec_normal, X_normal(end, :), 'linewidth', 1.2)
hold on
plot(tvec, X_ss(end-1,:),'--',tvec, X_ss(end,:), '--', 'linewidth', 1.2)
legend('vodd-normal','veven-normal','vodd-periodic','veven-periodic')
hold off;
ylabel('Voltage (V)')
xlabel('time (ns)')

