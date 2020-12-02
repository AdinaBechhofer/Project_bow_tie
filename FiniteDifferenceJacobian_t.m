function J = FiniteDifferenceJacobian_t(f,x,p,u,b,varargin)
epsilon = 0.01;
J = zeros(length(x), length(x));
for k=1:length(x)
    ek = zeros(length(x),1);
    ek(k) = 1;
    if ~isempty(varargin)
        t = varargin{1};
        J(:,k) = (f(x+epsilon*ek,p,u,b,t) - f(x,p,u,b,t))/epsilon;
    else 
        J(:,k) = (f(x+epsilon*ek,p,u,b) - f(x,p,u,b))/epsilon;
    end
end 
end 