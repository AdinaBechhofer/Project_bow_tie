function [X, t] = TrapMethodDynamicStep(x0,p,u_func,b,fJfhand, t_start, t_end, deltf_max)

X(:,1) = x0;
t(1) = t_start;
t_curr = t_start;
iter = 1;
while t_curr <t_end
    
    % Compute gamma
    %disp(size(u_func.func1(t_curr)))
    %disp(size(u_func.func2(t_curr)))
    u = u_func(t_curr);
    f0 = fJfhand(x0,p,u,b,t_curr); 
    dt = deltf_max/max(abs(f0));
    %f0 = fJfhand(x0); % our input is time independent
    gamma = x0 + (dt/2)*f0;
    % Netwon iteration
    itpause = 1;
    % initial guess for Newton is just x from prev timestep
    xt = newtonNd(@FJFTrap,x0,p,u_func(t_curr + dt),b,t_curr,gamma, dt,fJfhand);
    X(:,iter+1) = xt;   
    %diff = norm(X(:,l)-X(:,l+1))
    x0 = xt;
    t_curr = t_curr +dt;
    t(iter+1) = t_curr;
    iter = iter+1;
end