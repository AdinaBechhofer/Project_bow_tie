function [F JF] = ShootingMethod(x0,fhand,u)
epsilon = 0.01;
T=u.period;
tvec = 0:0.001:T;
Xtrap = TrapMethod(x0,fhand,tvec);
xt = Xtrap(end);
F = xt - x0;
JF = zeros(length(x0));
for k=1:length(x0)
    ek = zeros(length(x0),1);
    ek(k) = 1;
    xt2 = TrapMethod(x0+epsilon*ek,fhand, tvec);
    JF(:,k) = 1/epsilon * (xt2(end)-xt);
end 
JF = JF - eye(length(x0));
