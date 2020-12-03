function f= eval_f2(x,p,u)
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
taby = p.taby;
invC = p.invC;
CG = p.CG;
Ivec = zeros(2*num_bowties,1);
% Jn = u.jnano;
% Ivec(1:2:end) = u.vEmitter/Re - Jn;
% Ivec(2:2:end) = u.vCollector/Rc + Jn;
for i = 1:2*num_bowties
    if mod(i,2)==1
        Ivec(i) = u.vEmitter/Re - A*Jnano(phi,beta*(x(i)-x(i+1))/d,ROC,taby);
    else
        Ivec(i) = u.vCollector/Rc + A*Jnano(phi,beta*(x(i-1)-x(i))/d,ROC,taby);
        %disp(A*Jnano(phi,beta*(x(i-1)-x(i))/d,ROC,taby))
    end
    
end
f = CG*x+invC*Ivec;
disp(size(CG))
% % Stamp G here if non-linear
% f = zeros(2*num_bowties, 1);
% for i =1:2*num_bowties
%     if mod(i,2) ==1
%     f(i) = ((-x(i)/Re - Jn +u.vEmitter)/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i+1)/Rc + Jn +u.vCollector))/(1-Cec^2/(Cec+Cp)^2);
%     else
%     f(i) = ((-x(i)/Rc + Jn +u.vCollector)/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i-1)/Re - Jn +u.vEmitter))/(1-Cec^2/(Cec+Cp)^2);
%     end 
% end 