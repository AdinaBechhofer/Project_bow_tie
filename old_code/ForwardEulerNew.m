function X = ForwardEulerNew(fhand, x0,p,u,tvec)
X(:,1) = x0;
for n=1:length(tvec)-1
   dt = tvec(n+1) - tvec(n);
   f = fhand(X(:,n),tvec(n), p,u, p.b);
   X(:,n+1)= X(:,n) +  (dt * f);
end