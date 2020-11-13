clear;

p.gamma = 0.1;
p.km = 0.1;
p.ka=0.1;
p.N=500;
p.del_z = 1/(1+p.N);

A = zeros(p.N, p.N);
for i =1:p.N
    A(i,i) = -2*p.km/p.del_z^2 - p.ka;
    if i>1
        A(i, i-1) = p.km/p.del_z^2;
    end
    if i<p.N
        A(i, i+1) = p.km/p.del_z^2;
    end
end
A = 1/p.gamma*A;
b = zeros(p.N, 1);
b(1) = 1;
c = zeros(p.N, 1);
c(end) = 1;

[V,D] = eig(A);
b_tild = V'*b;
c_tild = V'*c;
%f= 1/p.gamma*A*x +b*u;
%y = c'*x;
search_factors = abs(b_tild.*c_tild)./abs(diag(D));
search_factors2 = abs(b_tild.*c_tild);
search_factors3 = abs(1./diag(D));
[out,idx] = sort(abs(search_factors3));

A_hat_2 = D(idx(end-1:end),idx(end-1:end));
b_hat_2 = b_tild(idx(end-1:end));
c_hat_2 = c_tild (idx(end-1:end));
V_q5 = V(:, idx(end-4:end));
A_hat_5 = D(idx(end-4:end),idx(end-4:end));
b_hat_5 = b_tild(idx(end-4:end));
c_hat_5 = c_tild (idx(end-4:end));
V_q10 = V(:, idx(end-9:end));
A_hat_10 = D(idx(end-9:end),idx(end-9:end));
b_hat_10 = b_tild(idx(end-9:end));
c_hat_10 = c_tild (idx(end-9:end));


sys_org = ss(A,b,c',0);
sys2 = ss(A_hat_2,b_hat_2,c_hat_2',0);
sys5 = ss(A_hat_5,b_hat_5,c_hat_5',0);
sys10 = ss(A_hat_10,b_hat_10,c_hat_10',0);
% Log-spaced frequency vector
freq = logspace(-3,2,100);

% Bode plot
bodePlot(freq,sys_org,sys2, sys5, sys10);

% Step response
figure
step(sys_org, sys2, sys5, sys10)
tic;
x0 = zeros(p.N, 1);
p.A = A;
p.b = b;
p.c = c;
p.func = @heat_bar;
t = 0:0.01:10;
u = ones(size(t));
X = trapezoidal2(x0, p, u, t);
Y = p.c'*X;
time_full = toc;
tic;
x02 = zeros(2, 1);
p2.A = A_hat_2;
p2.b = b_hat_2;
p2.c = c_hat_2;
p2.func = @heat_bar;
X2 = trapezoidal2(x02, p2, u, t);
Y2 = p2.c'*X2;
time_2 = toc;
tic;

x05 = zeros(5, 1);
p5.A = A_hat_5;
p5.b = b_hat_5;
p5.c = c_hat_5;
p5.func = @heat_bar;
X5 = trapezoidal2(x05, p5, u, t);
Y5 = p5.c'*X5;
time_5 = toc;
tic;

x010 = zeros(10, 1);
p10.A = A_hat_10;
p10.b = b_hat_10;
p10.c = c_hat_10;
p10.func = @heat_bar;
X10 = trapezoidal2(x010, p10, u, t);
Y10 = p10.c'*X10;
time_10 = toc;

figure;
plot(t, Y, t,Y2,'.', t, Y5,'.', t, Y10,'.');
legend('full', 'q=2', 'q=5', 'q=10')
xlabel('t')
ylabel('y')

figure;
plot(t, abs(Y2 - Y), t,abs(Y5 - Y), t,abs(Y10 - Y))
legend('Error q=2', 'q=5', 'q=10')
xlabel('t')
ylabel('|y(t) - y_q (t)|')

disp(['Time full :' num2str(time_full)])
disp(['Time q=2 :' num2str(time_2)])
disp(['Time q=5 :' num2str(time_5)])
disp(['Time q=10 :' num2str(time_10)])