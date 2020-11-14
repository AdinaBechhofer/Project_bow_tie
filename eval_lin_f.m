function f_approx = eval_lin_f(x,p,u)

[A, B] = linearize_in_x(@eval_f3, x0, p, b)
f_approx = A*x+B*[1; u.vEmitter; u.vCollector];