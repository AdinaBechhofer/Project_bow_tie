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
u =  [v1/p.REmitter;v2/p.RCollector];

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
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) + p.Ccoupling;
        end
    end
end
% figure;
% spy(C)
p.invC = inv(C);
p.CG = p.invC*G;

x0 = zeros(2*p.NumBowties, 1);
unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);
p.b = b;
%b(:, 1:2:10) = 1;
%b(:, 11:end) = 0;
c = eye(2*p.NumBowties);
x_ss = newtonNd(@fjbowtie,x0,p,u,b);
[A,B] = linearization(@eval_f3, x_ss, p, u,b);

[V,D] = eig(A);
b_tild = V'*b;
c_tild = V'*c;

freqs = [0.001, 0.003, 0.01, 0.03, 0.1,0.3, 1,3,10,30,100,300,1000,3000, 10000];
H_mag_1 = zeros(size(freqs));
H_phase_1 = zeros(size(freqs));
H_mag_2 = zeros(size(freqs));
H_phase_2 = zeros(size(freqs));


for q =1:length(freqs)
    s_minus_lambd = diag(1./(freqs(q) - diag(D)));
    H = c_tild*s_minus_lambd*b_tild;
    H_mag_1(q) = abs(H(1,1));
     H_phase_1(q) = angle(H(1,1));
     H_mag_2(q) = abs(H(2,1));
     H_phase_2(q) = angle(H(2,1));
end 

% freqs = [0.001, 0.003, 0.01, 0.03, 0.1,0.3, 1,3,10,30,100,300,1000,3000, 10000]
% H_mag_1 = zeros(size(freqs));
% H_phase_1 = zeros(size(freqs));
% H_mag_2 = zeros(size(freqs));
% H_phase_2 = zeros(size(freqs));
% for i =1:length(freqs)
%     omega = freqs(i);
%     H = c'*inv(1j*omega*eye(2*p.NumBowties) - A)*b;
%     H_mag_1(i) = abs(H(1,1));
%     H_phase_1(i) = angle(H(1,1));
%     H_mag_2(i) = abs(H(2,1));
%     H_phase_2(i) = angle(H(2,1));
% end 


figure;
subplot(2,1,1)
semilogx(freqs, mag2db(H_mag_1), freqs, mag2db(H_mag_2))
xlabel('Sine frequency')
ylabel('dB response')
title('Magnitude response')
legend('|H_{odd}|', '|H_{even}|')
subplot(2, 1, 2)
semilogx(freqs, H_phase_1, freqs, H_phase_2)
xlabel('Sine frequency')
ylabel('radians')
title('Phase response')
legend('<H_{odd}', '<H_{even}')