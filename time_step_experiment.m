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

x0 = zeros(2*p.NumBowties, 1);
unitb = [1, 0; 0, 1];
b = repmat(unitb,p.NumBowties,1);
p.b = b;

t_stop = 8;
t_start = 0;
del_ts = [0.04, 0.02, 0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001];
difference = zeros(length(del_ts), 80);
% stable = 0.12
del_t= 0.08;
tvec1 = t_start:del_t:t_stop;
% Implement forward euler
X1 = ForwardEulerNew(@eval_f4,x0,p,u,tvec1);

hold on 
for q=1:length(del_ts)
    tvec = t_start:del_ts(q):t_stop;
    X_curr = ForwardEulerNew(@eval_f4,x0,p,u,tvec);
    tsin = timeseries(X_curr',tvec);
    tsout = resample(tsin, tvec1);
    diff = abs(tsout.Data' - X1);
    plot(tvec1, diff(5,:), tvec1, diff(6, :), tvec1, diff(5+4*p.col,:), tvec1, diff(6+4*p.col,:))
    X1 = tsout.Data';
end 
X_ref = X_curr;
tvec_ref = tvec;
legendcell = strcat('\Delta t =',string(num2cell(repelem(del_ts, 4))));
legend(legendcell)
hold off

% max \del t before instability 
del_t = 0.12;
tvec1 = t_start:del_t:t_stop;
X1 = ForwardEulerNew(@eval_f4,x0,p,u,tvec1);
tsin_ref = timeseries(X_ref',tvec_ref);
tsout = resample(tsin_ref, tvec1);
diff = abs(tsout.Data' - X1);
figure;
plot(tvec1, diff(5,:), tvec1, diff(6, :), tvec1, diff(5+4*p.col,:), tvec1, diff(6+4*p.col,:));
% the error here is 0.5 of the final value of the simulation. Even though
% this isn't unstable, the error is wayy to high for us to accept.  

% delt t = 0.001 produces error of 6.4e-3 which is acceptable error. It's
% 0.13% error
tvec = t_start:0.001:t_stop;
X_curr = ForwardEulerNew(@eval_f4,x0,p,u,tvec);
tsout = resample(tsin_ref, tvec);
diff = abs(tsout.Data' - X_curr);
figure;
plot(tvec, diff(5,:), tvec, diff(6, :), tvec, diff(5+4*p.col,:), tvec, diff(6+4*p.col,:))
% 0.0005 << 0.12 (more than 2 orders of magnitude diff)
% max \del t before instability 

% let's check trap 
del_t = 0.02;
% del t = 0.02 results in err= 6.4e-3 which is 0.13% error
tvec1 = t_start:del_t:t_stop;
[X, tf_prod] = TrapMethod(x0,@(xp,tp)fjbowtie2(xp,tp,p,u,b), tvec1);
tsout = resample(tsin_ref, tvec1);
err = abs(tsout.Data' - X);
figure;
plot(tvec1, err(5,:), tvec1, err(6, :), tvec1, err(5+4*p.col,:), tvec1, err(6+4*p.col,:));
hold on 
plot(tvec1, X(5,:), tvec1, X(6, :), tvec1, X(5+4*p.col,:), tvec1, X(6+4*p.col,:), 'linewidth', 1.4);
title('Trap error and result')
hold off

disp(max(tf_prod))
max_dtf = max(tf_prod)/5;

[X_dynamic, t_points] = TrapMethodDynamicStep(x0,@(xp,tp)fjbowtie2(xp,tp,p,u,b),t_start, t_stop, max_dtf);
figure;
plot(t_points, X_dynamic(5,:), t_points, X_dynamic(6, :))
hold on 
tsout = resample(tsin_ref, t_points);
err = abs(tsout.Data' - X_dynamic);
plot(t_points, err(5,:), t_points, err(6, :))
plot(tvec_ref, X_ref(5,:), tvec_ref, X_ref(6, :))