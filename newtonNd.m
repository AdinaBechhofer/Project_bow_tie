function [x0,f] = newtonNd(fJfhand,x0,p, u, b, varargin)
% function newton1d(fhand,x0,itpause)
% 
% INPUTS:
%   fhand - function handle
%   Jhand - jacobian function handle    
%   x0    - initial guess
%   N,V   - input parameters for poisson
%   itpause - parameter for plotting
% 
% Use Newton's method to solve the nonlinear function
% defined by function handle fhand with initial guess x0.  
%
% itpause is parameter for plotting and defines the 
% number of Newton steps that are plotted sequentially
% pauses between sub-steps

if nargin<3
    error('Must provide three input arguments.  Type ''help newton1d'' for details');
end

tol=1e-10;          % convergence tolerance

% if isempty(varargin)
%     maxIters=500;       % max # of iterations
% else
%     maxIters = varargin{1}
% end
x00=x0;             % initial guess
maxIters=500;

% Newton loop
for iter=1:maxIters
    if ~isempty(varargin)
        [f, J] = fJfhand(x0, p, u, b, varargin);
    else 
        [f, J] = fJfhand(x0, p, u, b);
    end 
    %f = fhand(x0,N,V);              % evaluate function
    %J = Jhand(x0,N,V);              % evaluate jacobian
    dx=-J\f;                    % solve linear system
    nf(iter)=norm(f,inf);            % norm of f at step k+1
    ndx(iter)=norm(dx,inf);          % norm of dx at step k+1
    
    x(:,iter)=x0+dx;              % solution x at step k+1
    x0=x(:,iter);                 % set value for next guess

    if ndx(iter) < tol,          % check for convergence

        %fprintf('Converged to mean x=%4.12e in %d iterations\n',mean(x0),iter);
        break; 
    end
end

if iter==maxIters, % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter); 
end
