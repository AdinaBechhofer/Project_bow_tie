function f= eval_f(x, p, u)
num_bowties = p.NumBowties;
Re = p.REmitter;
Rc = p.RCollector;
A = p.Area;
beta = p.Beta;
d = p.Distance;
phi = p.workFunction;
Cec = p.CemitterCollector;
Cp = p.Cparasitic;
f = zeros(2*num_bowties, 1);
for i =1:2*num_bowties
    if mod(i,2) ==0
    f(i) = ((-x(i)/Re - A*jnano(beta*(x(i) - x(i+1))/d) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i+1)/Rc +A*jnano(beta*(x(i) - x(i+1))/d)+u(i+1)))/(1-Cec^2/(Cec+cp)^2);
    else
    f(i) = ((-x(i)/Rc + A*jnano(beta*(x(i-1) - x(i))/d) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i-1)/Re -A*jnano(beta*(x(i-1) - x(i))/d)+u(i-1)))/(1-Cec^2/(Cec+cp)^2);
    end 
end 