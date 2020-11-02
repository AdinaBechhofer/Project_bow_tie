function J = FiniteDifferenceJacobian(f, x, p, u)
epsilon = 0.01;
J = zeros(length(x), length(x));
for k=1:length(x)
    ek = zeros(length(x),1);
    ek(k) = 1;
    J(:,k) = (f(x+epsilon*ek, p,u) - f(x, p,u))/epsilon;
end 
end 