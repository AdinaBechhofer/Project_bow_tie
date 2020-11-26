%clear all;
close all;

p.NumBowties = 20;
% number of rows
p.row = 5;
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

v1 = 5;
v2 = 1;
u.vEmitter =  v1/p.REmitter;
u.vCollector = v2/p.RCollector;
u.ampE = 0.5;
u.ampC = 0.5;
u.phaseE = 0;
u.phaseC = 0;

unitb = [1, 0, 0, 0; 0, 1, 0, 0];
b = repmat(unitb,p.NumBowties,1);

% excite only the left most row, i.e. add sinusodial for those
for i = 1:2*p.row
    if mod(i,2) == 1
        b(i,3) = 1;
    else
        b(i,4) = 1;
    end
end

x0 = ones(2*p.NumBowties, 1);

t_stop = 1;
t_start = 0;
timestep = 0.001;

tvec = t_start:timestep:t_stop;
u.period = 20*timestep;

% %u.amp = 0;
% %X0 = ForwardEuler_t(@fjbowtie,x0,p,u,tvec,b);
% %X0 = TrapMethod(x0,@(xp,tp)fjbowtie(xp,tp,p,u,b), tvec);

% %X = ForwardEuler_t(@fjbowtie,x0,p,u,tvec,b);
X = TrapMethod(x0,@(xp,tp)fjbowtie(xp,tp,p,u,b), tvec);


figure(1);
plot(tvec,X(1,:),'--')
hold on;
plot(tvec,X(2,:))

x0_new = newtonS(@(xp,tp)fjbowtie(xp,tp,p,u,b),x0,u,1);
xT = TrapMethod(x0_new,@(xp,tp)fjbowtie(xp,tp,p,u,b),tvec);

plot(tvec, xT(1,:),'--',tvec, xT(2,:))
legend('vodd-before','veven-before','vodd-periodic','veven-periodic')

hold off;
% plot(tvec,X(6,:), tvec,xT(6,:))
% hold off
% ylabel('v(t)')
% xlabel('t')
% legend('vodd-before','vodd-steady-state','veven-before','veven-steady-state')


% plot(X(5+4*p.col,:))
% plot(X(6+4*p.col,:))
% plot(X(5+2*(p.row-1)*p.col, :))
% plot(X(6+2*(p.row-1)*p.col, :))
% legend('vodd1','veven1','vodd2','veven2','vodd3','veven3')
% title('voltages over time','FontSize', 8)
% hold off;




