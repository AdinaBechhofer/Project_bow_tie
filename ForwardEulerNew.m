function X = ForwardEulerNew(x0,p,u,t_start,t_stop,timestep, taby)
X(:,1) = x0;
t(1) = t_start;
for n=1:ceil((t_stop-t_start)/timestep)
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   %u = feval(eval_u, t(n));
   f = eval_f(X(:,n),p,u, taby);
   %f = eval_f_simSIR(X(:,n), p);
   X(:,n+1)= X(:,n) +  (dt * f);
end