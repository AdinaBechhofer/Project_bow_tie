function [X, t] = TrapMethodDynamicStep(x0, fJfhand, tstart, tend, deltf_max)

X(:,1) = x0;
t(1) = tstart;
t_curr = tstart;
iter = 1;
while t_curr <tend
    
    % Compute gamma
    f0 = fJfhand(x0,t_curr); 
    dt = deltf_max/max(f0);
    %f0 = fJfhand(x0); % our input is time independent
    gamma = x0 + (dt/2)*f0;
    % Netwon iteration
    itpause = 1;
    % initial guess for Newton is just x from prev timestep
    xt = newtonNd(fJfhand,gamma,x0,t_curr,dt,itpause);
    X(:,iter+1) = xt;   
    %diff = norm(X(:,l)-X(:,l+1))
    x0 = xt;
    t_curr = t_curr +dt;
    t(iter+1) = t_curr;
    iter = iter+1;
end