function [F, JF] = ShootingMethod(x0,p,U,b, fJfhand, tvec)
epsilon = 0.001;
%Xtrap = TrapMethod(x0,p,U,b, fJfhand, tvec);
Xeul = ForwardEuler_t(fJfhand,x0,p,U,b,tvec);
%xt = Xtrap(end);
xt = Xeul(:, end);
F = xt - x0;
JF = zeros(length(x0));
for k=1:length(x0)
    ek = zeros(length(x0),1);
    ek(k) = 1;
    %xt2 = TrapMethod(x0+epsilon*ek,p,U,b, fJfhand, tvec);
    xt2 = ForwardEuler_t(fJfhand,x0+epsilon*ek,p,U,b,tvec);
    JF(:,k) = 1/epsilon * (xt2(:, end)-xt);
end 
JF = JF - eye(length(x0));
