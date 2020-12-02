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
u = [v1/p.REmitter; v2/p.RCollector];

x0 = zeros(2*p.NumBowties, 1);
unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);
p.b = b;

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

t_stop = 8;
t_start = 0;
del_ts = [0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001, 0.00005, 0.000025, 0.00001];
difference = zeros(length(del_ts), 80);
% stable = 0.12
del_t= 0.02;
tvec1 = t_start:del_t:t_stop;
% Implement forward euler
X1 = ForwardEuler_t(@eval_f3,x0,p,u,b,tvec1);
comp_time = zeros(1, length(del_ts));
hold on 
for q=1:length(del_ts)
    tic;
    tvec = t_start:del_ts(q):t_stop;
    X_curr = ForwardEuler_t(@eval_f3,x0,p,u,b,tvec);
    comp_time(q) = toc;
    tsin = timeseries(X_curr',tvec);
    tsout = resample(tsin, tvec1);
    diff = abs(tsout.Data' - X1);
    %plot(tvec1, diff(5,:), tvec1, diff(6, :), tvec1, diff(5+4*p.col,:), tvec1, diff(6+4*p.col,:))
    plot(tvec1, max(diff, [], 1));
    X1 = tsout.Data';
end 
X_ref = X_curr;
tvec_ref = tvec;
legendcell = strcat('\Delta t =',string(num2cell(del_ts)));
title('Convergence of ForwardEuler with decreasing \Delta t')
ylabel('error')
xlabel('time')
legend(legendcell)
hold off

% max \del t before instability 
del_t = 0.12;
tvec1 = t_start:del_t:t_stop;
X1 = ForwardEuler_t(@eval_f3,x0,p,u,b,tvec1);
tsin_ref = timeseries(X_ref',tvec_ref);
tsout = resample(tsin_ref, tvec1);
diff = abs(tsout.Data' - X1);
figure;
plot(tvec1, max(diff, [], 1));
title('Error of ForwardEuler when ran with largets \Delta t for which it is stable')
ylabel('error')
xlabel('time')
% the error here is 0.5 of the final value of the simulation. Even though
% this isn't unstable, the error is wayy to high for us to accept.  

% delt t = 0.0007 produces error of 4.3e-3 which is acceptable error. It's
% 0.17% error
del_t_opt = 0.0005;
tic;
tvec_opt = t_start:del_t_opt:t_stop;
X_opt = ForwardEuler_t(@eval_f3,x0,p,u,b,tvec_opt);
tsout = resample(tsin_ref, tvec_opt);
diff = abs(tsout.Data' - X_opt);
figure;
plot(tvec_opt, max(diff, [], 1))
forward_euler_fixed_step_time = toc;
% 0.0005 << 0.12 (more than 2 orders of magnitude diff)
% max \del t before instability 

% let's check trap 
del_t = 0.02;
% del t = 0.02 results in err= 6.4e-3 which is 0.13% error
tic;
tvec1 = t_start:del_t:t_stop;
[X, tf_prod] = TrapMethod(x0,p,u,b,@fjbowtie, tvec1);
fixed_step_time = toc;
tsout = resample(tsin_ref, tvec1);
err = abs(tsout.Data' - X);
figure;
plot(tvec1, max(err, [], 1));
hold on 
plot(tvec1, X(5,:), tvec1, X(6, :), tvec1, X(5+4*p.col,:), tvec1, X(6+4*p.col,:), 'linewidth', 1.4)
title('Trap error and result')
hold off

disp(max(tf_prod))
max_dtf = max(tf_prod)/5;
tic;
[X_dynamic, t_points] = TrapMethodDynamicStep(x0,p,u,b,@fjbowtie,t_start,t_stop, max_dtf);
dynamic_time = toc;
figure;
plot(t_points, X_dynamic(5,:), t_points, X_dynamic(6, :), 'linewidth', 1.5)
hold on 
tsout = resample(tsin_ref, t_points);
err = abs(tsout.Data' - X_dynamic);
plot(t_points, err(5,:), t_points, err(6, :))
plot(tvec_ref, X_ref(5,:), tvec_ref, X_ref(6, :))
stem(t_points,  X_dynamic(5,:), 'r')
plot(tvec1, zeros(size(tvec1)), 'b.')
hold off 
title('Dynamic step size')
xlabel('time (fs)')
ylabel('Voltage (V)')