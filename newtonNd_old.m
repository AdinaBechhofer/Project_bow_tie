function x0 = newtonNd_old(fhand,x0,p,u,b,itpause,varargin)
% function newton1d(fhand,x0,itpause)
% 
% INPUTS:
%   fhand - function handle
%   x0    - initial guess
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

ftol=1e-10;          % convergence tolerance
dxtol = 1e-8;
maxIters=500;       % max # of iterations
x00=x0; % initial guess
p.x00 = x00;
X= zeros(length(x0),1);
% Newton loop
for iter=1:maxIters
    %disp(p)
    [f, J] = fhand(x0,p.t,p,u,b);    %x,t,p,u,b
    dx=J\(-f);                    % solve linear system
    nf(iter)=max(abs(f));            % norm of f at step k+1
    ndx(iter)=max(abs(dx));          % norm of dx at step k+1
    X(:,iter)=x0+dx;              % solution x at step k+1
    x0=x0+dx;                 % set value for next guess
    if nf(iter) <ftol && ndx(iter) < dxtol,          % check for convergence
        % check for convergence
        fprintf('Converged in %d iterations\n',iter);
        break; 
    end
end

if iter==maxIters, % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter); 
end

% stuff for plotting
 X=[x00,X];
% xmax=max(abs(x));
% xrange=linspace(-xmax,xmax,5000);
% [frange Jrange]=fhand(xrange, q)

% [fg Jg]=fhand(x, q);
 iters=1:iter;
figure(2)
subplot(1,2,1)
semilogy(iters, nf);
xlabel('iteration')
ylabel('||f||_{\infty}')
subplot(1,2,2)
semilogy(iters, ndx)
ylabel('||dx||_{\infty}')
xlabel('iterations')

figure(7)
subplot(1,2,1)
semilogy(iters, nf/nf(1));
xlabel('iteration')
ylabel('$||f(x^k)||_{\infty}/||f(x^0)||_{\infty}$', 'interpreter','latex')
title('Direct solve')
subplot(1,2,2)
semilogy(iters, ndx/ndx(1))
ylabel('$||dx^k||_{\infty}/||dx^0||_{\infty}$', 'interpreter','latex')
xlabel('iterations')
title('Direct solve')

% z_points = 0+p.del_z:p.del_z:1-p.del_z;
% figure(1);
% plot(z_points,x00)
% hold on
% for k =2:iter+1
%     plot(z_points,X(:,k))
% end
% legendCell = strcat('iteration ',string(num2cell(0:iter)));
%  legend(legendCell);
% hold off
% plot a few things
if ~isempty(varargin) && strcmp(varargin{1},'off')
else
%     plot_newton1d(x,fg,iters,ndx,xrange,frange,itpause);
end