function J = FiniteDifferenceJacobian_t(f, x, t, p,u,b)
epsilon = 0.01;
J = zeros(length(x), length(x));
for k=1:length(x)
    ek = zeros(length(x),1);
    ek(k) = 1;
    J(:,k) = (f(x+epsilon*ek, t,p,u,b) - f(x, t,p,u,b))/epsilon;
end 
end 