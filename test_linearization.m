clear all;
%close all;

% hint 1
x0 = 3;
%u0.vEmitter= 3;
%u0.vCollector = 0.25;
u = [3; 0.25];
xrange = -2:0.1:10;
[A,B] = linearization(@f_simple, x0, 0,0, u,0);
f_approx = A*xrange + B*[1; u];
figure;
plot(xrange, f_simple(xrange,0,0, u,0), xrange, f_approx)
hold on 
x0 = 6;
[A,B] = linearization(@f_simple, x0, 0,0, u, 0);
f_approx = A*xrange + B*[1; u];
plot(xrange, f_approx)
xlabel('x')
ylabel('f(x,p,u)')
title('Test tangent')
hold off 

% hint 2
x1 = zeros(100,1);
[A1,B1] = linearization(@f_linear, x1, 0,0, u, 0);
x2 = 5000*rand(100,1);
[A2,B2] = linearization(@f_linear, x2, 0,0, u, 0);

figure;
surf(A1-A2)
title('A1 - A2. Check independence from x')

%hint 3
x2 = zeros(100,1);
x2(1:3:end)= -5;
x2(2:3:end) = 1.7;
x2(2:5:end) = 3;
[A,B] = linearization(@f_special, x2, 0, u, 0);

clear all

p.NumBowties = 100;
p.REmitter = 1; % 1 ohm 
p.RCollector = 1; % 1 ohm  
p.Area = 50000e-18; % 5 nm^2 
p.Beta = 25; % enhancement factor 
p.Distance = 10; % 10 nm
p.workFunction = 5.1; % work function of gold 
p.CemitterCollector = 2; % nano farad
p.Cparasitic = 0.05; % 0.1 nano farad
p.Radius = 10; % 1 nm
p.taby = csvread('rspa20140811supp3.csv');    
% number of rows
p.row = 10;
% number of columns 
p.col = 10;
p.Ccoupling = 0.03;

v1 = 5;
v2 = 0;
u =  [v1/p.REmitter; v2/p.RCollector];

C = zeros(2*p.NumBowties);
G = zeros(2*p.NumBowties);
% Stamp C and invert; stamp G
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
        C(i,i-1)= C(i,i-1) -p.CemitterCollector;
        G(i,i) = -1./(p.RCollector);
        if i<= 2*p.row*(p.col-1)
            C(i,i) = C(i,i) + p.Ccoupling;
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) - p.Ccoupling;
        end
    end
end
p.invC = inv(C);
p.CG = p.invC*G;

x0 = zeros(2*p.NumBowties, 1);
t_stop = 15;
t_start = 0;
timestep = 0.01;

% Implement forward euler
X(:,1) = x0;
X_lin(:,1) = x0;
t(1) = t_start;
for n=1:ceil((t_stop-t_start)/timestep)
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   f = eval_f2(X(:,n),p,u);
   X(:,n+1)= X(:,n) +  (dt * f);
   [A,B] = linearization(@eval_f2, X_lin(:,n), p, u);
   f_approx = A*X_lin(:,n)+B*[1; u];
   X_lin(:,n+1)= X_lin(:,n) +  (dt * f_approx);
end 

figure;
plot(t, X(5,:), t, X_lin(5,:), 'k--')
hold on;
plot(t, X(6,:), t, X_lin(5,:), 'k--')
%plot(X(5+4*p.col,:))
%plot(X(6+4*p.col,:))
plot(t, X(5+2*(p.row-1)*p.col, :), t, X_lin(5+2*(p.row-1)*p.col, :), 'k--')
plot(t, X(6+2*(p.row-1)*p.col, :), t, X_lin(6+2*(p.row-1)*p.col, :), 'k--')
legend('vodd1','veven1','vodd_lin1','veven_lin1','vodd2','veven2','vodd_lin2','veven_lin2')
title('voltages over time','FontSize', 8)
hold off;

