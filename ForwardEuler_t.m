function [X, varargout] = ForwardEuler_t(fhand, x0,p,U,b, varargin)
X(:,1) = x0;

if length(varargin) ==1
    tvec = varargin{1};
    for n = 1:length(tvec)-1
       u = U(:, n);
       dt = tvec(n+1)-tvec(n);
       f = fhand(X(:,n),p,u,b, tvec(n));
       X(:,n+1)= X(:,n) +  (dt * f);
    end 
elseif length(varargin) >1
  if ~isstring(varargin{1})
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
  end
end
if nargout>1
    varargout{1} = tvec;
end
end 