function X = ForwardEuler_t(fhand, x0,p,u,tvec,b)
X(:,1) = x0;
t = tvec;

for n = 1:length(tvec)-1
   dt = tvec(n+1)-t(n);
   %u = feval(eval_u, t(n));
   f = fhand(X(:,n),t(n),p,u,b);
   %f = eval_f_simSIR(X(:,n), p);
   X(:,n+1)= X(:,n) +  (dt * f);
end