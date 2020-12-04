clear all;
close all;

p.NumBowties = 4;
% p.NumBowties = 25;
% number of rows
p.row = 2;
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
p.Cparasitic = 0.1; % 0.1 nano farad
% p.Radius = 1; % 1 nm?
p.Radius = 10; % 1 nm
p.taby = csvread('rspa20140811supp3.csv');    

% p.Ccoupling = 0.03;
p.Ccoupling = 0.08; %increase coupling
p.CcouplingV = 0.08; %increase coupling

p.jnano = 0; % to turn on/off non-linear jnano

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
amplitude = 0.5;
% amplitude = 10; % amplify effect of incident current
u = [p.v1; p.v2];
unitb = [1, 0, 0, 0; 0, 1, 0, 0];
b = repmat(unitb,p.NumBowties,1); % bias always on

% % turn on centre bowtie
% b(p.NumBowties, 3) = 1;
% b(round(p.NumBowties)+1, 4) = 1;

% turn on all bowties
b(1:2:end, 3) = 1;
b(2:2:end, 4) = 1;

x0 = zeros(2*p.NumBowties, 1);
tic;
% t_stop = 40;
t_stop_dc = 20;
t_start = 0;
timestep = 0.05;
tvec_normal = t_start:timestep:t_stop_dc;
tc = 25; % centre of pulse
fwhm = 5; % in ns
sigma = fwhm/2.35;
gauss = exp(-(tvec_normal-tc).^2/(sqrt(2*pi)*sigma));
U = [repmat(u,1,length(tvec_normal)); 
    amplitude*cos(2*pi*tvec_normal/period).* gauss; 
    amplitude/2*cos(2*pi*tvec_normal/period +0.03).* gauss];

% 
% U = [repmat(u,1,length(tvec_normal)); 
%     ones(1,length(tvec_normal)); 
%     ones(1,length(tvec_normal))];


%% linearize once near x0 = dc steady state
x_dc = newtonNd(@fjbowtie,x0,p,[u; 0; 0],b);

x0 = x_dc;

% solve Forward Euler from x_dc
X = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec_normal);

time_euler = toc

[VL, S, VR] = svd(X);
figure(1)
semilogx(diag(S))
ylabel('Singular values')
xlabel('index')
q = 5;
Vq = VL(:,1:q);

c1 = zeros(2*p.NumBowties,1);
c1(5) = 1; %consider 5th node
c2 = zeros(2*p.NumBowties,1);
c2(6) = 1;

c1hat = Vq'*c1;
c2hat = Vq'*c2;

xq0 = zeros(q, 1);

t_stop = 1; 
t_start = 0;
timestep = 0.05; 
X_lin_hat = zeros(length(xq0),ceil((t_stop-t_start)/timestep));


[A, B] = linearize_in_x(@eval_f3,x_dc,p,[u;1;1],b,tvec_normal(1));

Ahat = Vq'*A*Vq;
Bhat = Vq'*B;
Xo(:,1) = x_dc;
X_lin(:,1) = x_dc;
X_lin_hat(:,1) = Vq'*x_dc;
t(1) = t_start;

for n=1:ceil((t_stop-t_start)/timestep)
   u = U(:, n);
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   f = eval_f3(X(:,n),p,u,b, t(n));
   
   Xo(:,n+1)= Xo(:,n) +  (dt * f);
   f_approx = A*X_lin(:,n)+B*[1; u];
  
   X_lin(:,n+1)= X_lin(:,n) +  (dt * f_approx);
   f_approx_hat = Ahat*X_lin_hat(:,n)+Bhat*[1; u];
   X_lin_hat(:,n+1)= X_lin_hat(:,n) +  (dt * f_approx_hat);
end 

y1 = c1'*Xo;
y1_lin = c1'*X_lin;
y1_lin_hat5 = c1hat'*X_lin_hat;

y2 = c2'*Xo;
y2_lin = c2'*X_lin;
y2_lin_hat5 = c2hat'*X_lin_hat;

figure(2);
plot(t, y1, 'k.-')
hold on;
plot(t, y1_lin, 'b.-')
%plot(t, y1_lin_hat1,'r--')
%plot(t, y1_lin_hat2,'g--')
plot(t, y1_lin_hat5,'m--')
plot(t, y2, 'k.-')
plot(t, y2_lin, 'b.-')
%plot(t, y2_lin_hat1,'r--')
%plot(t, y2_lin_hat2,'g--')
plot(t, y2_lin_hat5,'m--')
legend('v1','v1-lin','v1-lin-q=5','v2','v2-lin','v2-lin-q=5')
title('voltages over time','FontSize', 8)
hold off;
ylabel('Voltage (V)')
xlabel('time (ns)')

% %% Linearize at every step
% % Implement forward euler with linearize system
% Xo(:,1) = x0;
% X_lin(:,1) = x0;
% X_lin_hat(:,1) = xq0;
% t(1) = t_start;
% for n=1:ceil((t_stop-t_start)/timestep)
%    dt = min(timestep, (t_stop-t(n)));
%    t(n+1)= t(n) + dt;
%    f = eval_f3(X(:,n),t(n),p,u,b);
%    Xo(:,n+1)= Xo(:,n) +  (dt * f);
%    [A, B] = linearize_in_x(@eval_f3,x0,t(n),p,u,b);
%    Ahat = Vq'*A*Vq;
%    Bhat = Vq'*B;
%    f_approx = A*X_lin(:,n)+B*[1; u.vEmitter; u.vCollector];
%    X_lin(:,n+1)= X_lin(:,n) +  (dt * f_approx);
%    f_approx_hat = Ahat*X_lin_hat(:,n)+Bhat*[1; u.vEmitter; u.vCollector];
%    X_lin_hat(:,n+1)= X_lin_hat(:,n) +  (dt * f_approx_hat);
% end 
% 
% y1 = c1'*Xo;
% y1_lin = c1'*X_lin;
% y1_lin_hat5 = c1hat'*X_lin_hat;
% 
% y2 = c2'*Xo;
% y2_lin = c2'*X_lin;
% y2_lin_hat5 = c2hat'*X_lin_hat;
% 
% figure;
% plot(t, y1, t, y1_lin)
% hold on;
% plot(t, y2,t, y2_lin)
% legend('vodd1','veven1','vodd_lin1','veven_lin1','vodd2','veven2','vodd_lin2','veven_lin2')
% title('voltages over time','FontSize', 8)
% hold off;
% 
% figure;
% plot(t, y1, 'k.-', t, y1_lin, 'b',t, y1_lin_hat,'r--')
% hold on;
% plot(t, y2, 'k.-', t, y2_lin, 'b',t, y2_lin_hat,'r--')
% 
% legend('v1','v1-lin','v1-lin-q=1','v2','v2-lin','v2-lin-q=1')
% title('voltages over time','FontSize', 8)
% hold off;
% 
% 
% figure;
% plot(t, y1, 'k.-', t, y1_lin, 'b',t, y1_lin_hat10,'r--',t, y1_lin_hat5,'c--', t, y1_lin_hat2, 'm--')
% hold on;
% plot(t, y2, 'k.-', t, y2_lin, 'b',t, y2_lin_hat10,'r--',t, y2_lin_hat5,'c--', t, y2_lin_hat2, 'm--')
% 
% legend('vodd','vodd-lin','vodd-lin-reduced','veven','veven-lin','veven-lin-reduced')
% title('voltages over time','FontSize', 8)
% hold off;