clear all;
close all;

p.NumBowties = 500;
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
%u.jnano = 1;

v1 = 5;
v2 = 0;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;

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

t_stop = 0.5;
t_start = 0;
timestep = 0.005;

% Generate solutions in time
tvec = t_start:timestep:t_stop;

unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);
% Implement forward euler
X = ForwardEuler_t(@eval_f3,x0,p,u,t_start,t_stop,timestep,b);
% 
% % plot some nodes
% figure(2);
% plot(X(5,:))
% hold on;
% plot(X(6,:))
% plot(X(5+4*p.col,:))
% plot(X(6+4*p.col,:))
% plot(X(5+2*(p.row-1)*p.col, :))
% plot(X(6+2*(p.row-1)*p.col, :))
% legend('vodd1','veven1','vodd2','veven2','vodd3','veven3')
% title('voltages over time','FontSize', 8)
% hold off;

%% Eigenmode truncation
c1 = zeros(2*p.NumBowties,1);
c1(5) = 1; %consider 5th node
c2 = zeros(2*p.NumBowties,1);
c2(6) = 1;

% consider times close to t=0 and linearize once at x0
[A, B] = linearize_in_x(@eval_f3,x0,0,p,u,b);
[V,D] = eig(A);
btilde = V\B;
c1tilde = V\c1;
c2tilde = V\c2;

Dinv = 1./diag(D);
residue = real(Dinv).*(abs(btilde.*c1tilde));
%residue = Dinv;
[sortedres, idx] = sort(residue,'ascend');

q = 1000;
Vq = V(:,idx(1:q));
Ahat = Vq'*A*Vq;
Bhat = Vq'*B;
c1hat = Vq'*c1;
c2hat = Vq'*c2;

xq0 = zeros(q, 1);

t_stop = 0.01; % just do 10 steps -- too slow!
t_start = 0;
timestep = 0.001; 
% Implement forward euler with linearize system
Xo(:,1) = x0;
X_lin(:,1) = x0;
X_lin_hat(:,1) = xq0;
t(1) = t_start;
for n=1:ceil((t_stop-t_start)/timestep)
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   f = eval_f3(X(:,n),t(n),p,u,b);
   Xo(:,n+1)= Xo(:,n) +  (dt * f);
   f_approx = A*X_lin(:,n)+B*[1; u.vEmitter; u.vCollector];
   X_lin(:,n+1)= X_lin(:,n) +  (dt * f_approx);
   f_approx_hat = Ahat*X_lin_hat(:,n)+Bhat*[1; u.vEmitter; u.vCollector];
   X_lin_hat(:,n+1)= X_lin_hat(:,n) +  (dt * f_approx_hat);
end 

y1 = c1'*Xo;
y1_lin = c1'*X_lin;
y1_lin_hat = c1hat'*X_lin_hat;

y2 = c2'*Xo;
y2_lin = c2'*X_lin;
y2_lin_hat = c2hat'*X_lin_hat;

% figure;
% plot(t, y1, t, y1_lin)
% hold on;
% plot(t, y2,t, y2_lin)
% legend('vodd1','veven1','vodd_lin1','veven_lin1','vodd2','veven2','vodd_lin2','veven_lin2')
% title('voltages over time','FontSize', 8)
% hold off;

figure;
plot(t, y1, 'k', t, y1_lin, 'b',t, y1_lin_hat,'r--')
hold on;
plot(t, y2, 'k', t, y2_lin, 'b',t, y2_lin_hat,'r--')
legend('vodd','vodd-lin','vodd-lin-reduced','veven','veven-lin','veven-lin-reduced')
title('voltages over time','FontSize', 8)
hold off;



