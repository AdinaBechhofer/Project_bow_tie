clear all;
close all;

p.NumBowties = 21*21;
% p.NumBowties = 49;
% number of rows
p.row = 21;
% p.row = 7;
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

p.jnano = 1; % to turn on/off non-linear jnano

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

% turn on centre bowtie
b(p.NumBowties, 3) = 1;
b(round(p.NumBowties)+1, 4) = 1;

% % turn on all bowties
% b(1:2:end, 3) = 1;
% b(2:2:end, 4) = 1;

t_stop = 20;
t_start = 0;
timestep = 0.05;
tvec = t_start:timestep:t_stop;

% Define pulse
tc = 10; % centre of pulse
fwhm = 5; % in ns
sigma = fwhm/2.35;
gauss = exp(-(tvec-tc).^2/(sqrt(2*pi)*sigma));
U = [repmat(u,1,length(tvec));
    amplitude*cos(2*pi*tvec/period).* gauss;
    amplitude/2*cos(2*pi*tvec/period +0.03).* gauss];

% Use Newton to find DC state
x0 = zeros(2*p.NumBowties, 1);
x_dc = newtonNd(@fjbowtie,x0,p,[u; 0; 0],b);

x0 = x_dc; % this is our new initial state

% Generate trajectories
X = [];
% b = repmat(unitb,p.NumBowties,1); % bias always on
% b(p.NumBowties, 3) = 1;
% b(p.NumBowties+1, 4) = 1;
% X1 = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec);
% X = [X, X1];
% X1 = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec);
% X = [X, X1];
% for n = 1:randsample(1:p.NumBowties,80)
for n = 1:p.NumBowties/3
    tvec = 0:0.05:10;
    b = repmat(unitb,p.NumBowties,1); % bias always on
    b(2*n-1, 3) = 1;
    b(2*n, 4) = 1;
    X1 = ForwardEuler_t(@eval_f3,x0,p,U,b,tvec);
    X1 = X1(:,1:2:end);
    X = [X, X1];
end

% Perform SVD
[VL, S, VR] = svd(X);
figure(1)
semilogy(diag(S),'LineWidth',2)
% ylim([1e-5,1e4])
% xlim([1,10])
ylabel('Singular values')
xlabel('index')

% Pick top q singular values
q = 450;
Vq = VL(:,1:q);
msetotal = 0;

%% Experimenting different pulses /inputs
% Define new pulse, random location
% testloc = 2*randsample(1:p.NumBowties,1) % this is even
% b = repmat(unitb,p.NumBowties,1); % bias always on
% b(testloc-1, 3) = 1;
% b(testloc, 4) = 1;
% testloc = 2*randsample(1:p.NumBowties,1) % this is even
b = repmat(unitb,p.NumBowties,1); % bias always on
b(p.NumBowties, 3) = 1;
b(p.NumBowties+1, 4) = 1;

t_stop = 20;
t_start = 0;
timestep = 0.005;
tvec = t_start:timestep:t_stop;
tc = 10; % centre of pulse
fwhm = 5; % in ns
period =1.5;
amplitude = 0.5;
u = [5; 0];
sigma = fwhm/2.35;
gauss = exp(-(tvec-tc).^2/(sqrt(2*pi)*sigma));
U = [repmat(u,1,length(tvec));
    amplitude*cos(2*pi*tvec/period).* gauss;
    amplitude/2*cos(2*pi*tvec/period +0.03).* gauss];

% track output for every single bowtie
% for ii = 1:2*p.NumBowties
%     c = zeros(2*p.NumBowties,1);
%     c(ii) = 1;
%     chat = Vq'*c;
%     X_lin_hat = zeros(q,ceil((t_stop-t_start)/timestep));
% 
%     % if linearize only x0, u and t are dummies (not used at all)
%     [A, B] = linearization(@eval_f3,x_dc,p,[u;1;1],b,tvec(1),'onlyx0');
% 
%     Ahat = Vq'*A*Vq;
%     Bhat = Vq'*B;
%     
%     Xo = ForwardEuler_t(@eval_f3,x_dc,p,U,b,tvec);
%     X_lin_hat = ForwardEuler_t(@eval_f3,Vq'*x_dc,p,U,b,'linearmodel',tvec,Ahat,Bhat);
%     
%     y = c'*Xo;
%     yhat = chat'*X_lin_hat;
%     MSE = sum((y-yhat).^2);
%     MSE = MSE/length(y)
%     msetotal = msetotal + MSE;
%     if ii == p.NumBowties
%         centreerror = MSE
%         y1_ori = y;
%         y1_red = yhat;
%     end
%     if ii == p.NumBowties + 1
%         centreerror2 = MSE
%         y2_ori = y;
%         y2_red = yhat;
%     end
% end
% msetotal/(2*p.NumBowties)
% centreerror
% centreerror2

% figure;
% subplot(2,1,1)
% plot(tvec, y1_ori)
% hold on
% plot(tvec, y1_red, '.')
% hold off
% legend('y1','y1-red')
% 
% subplot(2,1,2)
% plot(tvec, y2_ori)
% hold on
% plot(tvec, y2_red, '.')
% hold off
% legend('y2','y2-red')

% % consider random nodes
c1 = zeros(2*p.NumBowties,1);
% c1(5) = 1; %consider 5th node
c1(p.NumBowties) = 1; %consider 5th node
c2 = zeros(2*p.NumBowties,1);
% c2(6) = 1;
c2(p.NumBowties+1) = 1;
% 
c1hat = Vq'*c1;
c2hat = Vq'*c2;
% 
% X_lin_hat = zeros(q,ceil((t_stop-t_start)/timestep));
% 
% % if linearize only x0, u and t are dummies (not used at all)
[A, B] = linearization(@eval_f3,x_dc,p,[u;1;1],b,tvec(1),'onlyx0');

Ahat = Vq'*A*Vq;
Bhat = Vq'*B;

tic
Xo = ForwardEuler_t(@eval_f3,x_dc,p,U,b,tvec);
time_euler = toc

tic
X_lin = ForwardEuler_t(@eval_f3,x_dc,p,U,b,'linearmodel',tvec,A,B);
time_euler_lin = toc

tic
X_lin_hat = ForwardEuler_t(@eval_f3,Vq'*x_dc,p,U,b,'linearmodel',tvec,Ahat,Bhat);
time_euler_lin_red = toc

y1 = c1'*Xo;
y1_lin = c1'*X_lin;
y1_lin_hat = c1hat'*X_lin_hat;

y2 = c2'*Xo;
y2_lin = c2'*X_lin;
y2_lin_hat = c2hat'*X_lin_hat;

figure(2);
subplot(2,1,1)
plot(tvec, y1, 'k-','LineWidth',1)
hold on
plot(tvec, y1_lin, 'b.','MarkerSize',1)
% plot(tvec, y1_lin_hat1,'g.')
plot(tvec, y1_lin_hat,'m.','MarkerSize',1)
legend('v1','v1-lin','v1-lin-q=1','v1-lin-q=3')
title('voltages over time','FontSize', 8)
ylabel('Voltage (V)')
xlim([3,15])
hold off;
% 
subplot(2,1,2)
plot(tvec, y2, 'k-','LineWidth',1)
hold on
plot(tvec, y2_lin, 'b.','MarkerSize',1)
% plot(tvec, y2_lin_hat1,'g.')
plot(tvec, y2_lin_hat,'m.','MarkerSize',1)
xlim([3,15])
legend('v2','v2-lin','v2-lin-q=1','v2-lin-q=3')
hold off;
ylabel('Voltage (V)')
xlabel('time (ns)')

sum((y1_lin_hat - y1).^2)/length(y1)
sum((y2_lin_hat - y2).^2)/length(y2)