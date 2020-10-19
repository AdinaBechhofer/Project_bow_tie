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
ROC = p.Radius;
f = zeros(2*num_bowties, 1);
for i =1:2*num_bowties
    if mod(i,2) ==0
    f(i) = ((-x(i)/Re - A*jnano(phi,beta*(x(i) - x(i+1))/d, ROC) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i+1)/Rc +A*jnano(phi,beta*(x(i) - x(i+1))/d, ROC)+u(i+1)))/(1-Cec^2/(Cec+cp)^2);
    else
    f(i) = ((-x(i)/Rc + A*jnano(phi,beta*(x(i-1) - x(i))/d, ROC) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i-1)/Re -A*jnano(phi,beta*(x(i-1) - x(i))/d, ROC)+u(i-1)))/(1-Cec^2/(Cec+cp)^2);
    end 
end 