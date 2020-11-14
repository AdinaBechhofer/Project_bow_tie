function X = TrapMethod(x0, fJfhand, tvec)

X(:,1) = x0;
t = tvec;

for l = 1:length(tvec)-1
    
    dt = tvec(l+1)-tvec(l);
    
    % Compute gamma
    f0 = fJfhand(x0,t(l)); 
    %f0 = fJfhand(x0); % our input is time independent
    gamma = x0 + (dt/2)*f0;
    
    % Netwon iteration
    itpause = 1;
    % initial guess for Newton is just x from prev timestep
    xt = newtonNd(fJfhand,gamma,x0,t(l),dt,itpause);
    X(:,l+1) = xt;   
    %diff = norm(X(:,l)-X(:,l+1))
    
    x0 = xt;
    
  
end