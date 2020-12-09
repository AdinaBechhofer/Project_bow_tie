function test_demo(pulsestr)
% clear all;
close all;

p.NumBowties = 49;
% p.NumBowties = 25;
% number of rows
p.row = 7;
% p.row = 5;
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
% p.Cparasitic = 0.05; % 0.1 nano farad
p.Cparasitic = 0.3; % 0.1 nano farad
% p.Radius = 1; % 1 nm?
p.Radius = 10; % 1 nm
p.taby = csvread('rspa20140811supp3.csv');    

% p.Ccoupling = 0.03;
p.Ccoupling = 0.28; %increase coupling
p.CcouplingV = 0.28; %increase coupling
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
            C(i,i) = C(i,i) + p.Ccoupling ;
            C(i,i-2*p.row+1) = C(i,i-2*p.row+1) - p.Ccoupling;
        end 
        if mod(i-1,p.row*2)~=0 && mod(i+1,p.row*2)
            C(i,i) = C(i,i) + 2*p.CcouplingV ;
            C(i,i-2) = C(i,i-2) -p.CcouplingV;
            C(i,i+2) = C(i,i+2) -p.CcouplingV;
        end
    else
        C(i,i-1)= -p.CemitterCollector;
        G(i,i) = -1./(p.RCollector);
        if i<= 2*p.row*(p.col-1)
            G(i,i) = G(i,i) + p.Ccoupling;
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) - p.Ccoupling;
        end
        if mod(i-2,p.row*2)~=0 && mod(i,p.row*2)~=0
            C(i,i) = C(i,i) + 2*p.CcouplingV ;
            C(i,i-2) = C(i,i-2) -p.CcouplingV;
            C(i,i+2) = C(i,i+2) -p.CcouplingV;
        end
    end
end
% figure;
% spy(C)
p.invC = inv(C);
p.CG = p.invC*G;

p.v1 = 5;
p.v2 = 0;
period =1.5;
%amplitude = 0.5;
amplitude = 10; % amplify effect of incident current
u = [p.v1; p.v2];

[b, demopt] = load_pulse_position(p,pulsestr);

% [b, demopt] = load_pulse_position(p,'centre-only');
% [b, demopt] = load_pulse_position(p,'centre-neighbours');
% [b, demopt] = load_pulse_position(p,'left-side');
% [b, demopt] = load_pulse_position(p,'left-corner');

x0 = zeros(2*p.NumBowties, 1);
x_dc = newtonNd(@fjbowtie,x0,p,[u; 0; 0],b);
x0 = x_dc; % this is our new initial state

t_stop = 40;
t_start = 0;
timestep = 0.05;
tvec_normal = t_start:timestep:t_stop;
tc = 15; % centre of pulse
%tc = 25;
fwhm = 40; % in ns
%fwhm = 5; % in ns
sigma = fwhm/2.35;
gauss = exp(-(tvec_normal-tc).^2/(sqrt(2*pi)*sigma));
U = [repmat(u,1,length(tvec_normal)); 
    amplitude*cos(2*pi*tvec_normal/period).* gauss; 
    amplitude/2*cos(2*pi*tvec_normal/period +0.03).* gauss];

tic;
X_normal = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec_normal);
time_euler = toc;

tic
for n = demopt
    VisualizeNetwork(X_normal(:,n),p,U(:,n),pulsestr)
    drawnow
end
time_visual = toc;

% x0 = zeros(2*p.NumBowties, 1);
% tic;
% X_normal = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec_normal);
% time_euler = toc

% figure;
% subplot(2, 1, 1)
% plot(tvec_normal, X_normal(p.NumBowties, :), tvec_normal, X_normal(p.NumBowties+1, :), 'linewidth', 1.2)
% legend('centre-v1','centre-v2')
% ylabel('Voltage (V)')
% xlabel('time (ns)')
% 
% subplot(2, 1, 2)
% plot(tvec_normal, X_normal(p.NumBowties+4, :), tvec_normal, X_normal(p.NumBowties+5, :), 'linewidth', 1.2)
% legend('next nearest neighour-v1','next nearest neighbour-v2')
% ylabel('Voltage (V)')
% xlabel('time (ns)')


