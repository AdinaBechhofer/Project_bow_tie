p.R = 1;
p.C = 0.1;

x0 = 0;
u0 = 5;


t_stop = 5;
t_start = 0;
timestep = 0.01;

% Implement forward euler
X = ForwardEulerNewRC(x0,p,u0,t_start,t_stop,timestep);

figure;
plot(X)

legend('v')
title('voltages over time','FontSize', 8)
hold off;

