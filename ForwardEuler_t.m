function [X, varargout] = ForwardEuler_t(fhand, x0,p,U,b, varargin)
X(:,1) = x0;
if length(varargin) ==1
    tvec = varargin{1};
    tf_prod = zeros(1, length(tvec));
    for n = 1:length(tvec)-1
       u = U(:, n);
       dt = tvec(n+1)-tvec(n);
       f = fhand(X(:,n),p,u,b, tvec(n));
       tf_prod = max(dt*f);
       X(:,n+1)= X(:,n) +  (dt * f);
    end 
elseif length(varargin) >1
  if varargin{1}~='dynamic'
    t_start = varargin{1};
    t_stop = varargin{2};
    del_t = varargin{3};
    tvec(1) = t_start;
    for n=1:ceil((t_stop-t_start)/del_t)
       dt = min(del_t, (t_stop-tvec(n)));
       u = U(:, n);
       tvec(n+1)= tvec(n) + dt;
       f = fhand(X(:,n),p,u, b,tvec(n));
       X(:,n+1)= X(:,n) +  (dt * f);
    end
  else 
      % dyhnamic
    t_start = varargin{2};
    t_stop = varargin{3};
    max_prod = varargin{4};
    tvec(1) = t_start;
    t_curr = t_start;
    n = 1;
    while t_curr <t_stop
       f = fhand(X(:,n),p,U(t_curr), b,t_curr);
       dt = max_prod/max(abs(f));
       tvec(n+1) = tvec(n) +dt;
       X(:,n+1)= X(:,n) +  (dt * f);
       t_curr = t_curr + dt;
       n = n+1;
    end
      
  end
end
if nargout==2
    varargout{1} = tvec;
elseif nargout==3
    varargout{1} = tvec;
    varargout{2} = max(tf_prod);
end
end 