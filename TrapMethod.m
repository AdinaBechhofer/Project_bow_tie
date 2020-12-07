function [X, varargout] = TrapMethod(x0,p,U,b, fJfhand, tvec)

X(:,1) = x0;
t = tvec;
tf_prod = zeros(size(tvec));
for l = 1:length(tvec)-1
    
    dt = tvec(l+1)-tvec(l);
    % Compute gamma
    u = U(:, l);
    f0 = fJfhand(x0,p,u,b,t(l)); 
    %f0 = fJfhand(x0); % our input is time independent
    gamma = x0 + (dt/2)*f0;
    tf_prod(l) = max(abs(dt*f0));
    % Netwon iteration
    itpause = 1;
    % initial guess for Newton is just x from prev timestep
    xt = newtonNd(@FJFTrap,x0,p,U(:, l+1),b,t(l),gamma, dt,fJfhand);
    X(:,l+1) = xt;   
    %diff = norm(X(:,l)-X(:,l+1))
    
    x0 = xt;
end 
    if nargout>1
        varargout{1} = max(tf_prod);
    end 
  
end