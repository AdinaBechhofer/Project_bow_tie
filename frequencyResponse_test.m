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
p.t = 0;
v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;

x0 = zeros(2*p.NumBowties, 1);
unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);
p.b = b;
c = eye(2*p.NumBowties);
x_ss = newtonNd_old(@fjbowtie2,x0,p,u,b,8);
[A,B] = linearization(@eval_f4, x_ss,p.t, p, u,b);

freqs = [1,3,10,30,100,300,1000,3000, 10000]
H_mag = zeros(size(freqs));
H_phase = zeros(size(freqs));
for i =1:length(freqs)
    omega = freqs(i);
    H = c'*inv(1j*omega*eye(2*p.NumBowties) - A)*b;
    H_mag(i) = abs(H(2,1));
    H_phase(i) = angle(H(2,1));
end 

figure;
subplot(2,1,1)
semilogx(freqs, mag2db(H_mag))
xlabel('Sine frequency')
ylabel('dB response')
title('Magnitude response')
subplot(2, 1, 2)
semilogx(freqs, H_phase)
xlabel('Sine frequency')
ylabel('radians')
title('Phase response')