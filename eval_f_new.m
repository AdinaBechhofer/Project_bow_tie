function f= eval_f_new(x,p,u)
num_bowties = p.NumBowties;
num_in_col = p.Col;
num_in_row = p.Row;
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
Cwire = p.Cwire;
Rwire = p.Rwire;
Jn = u.jnano;

f = zeros(4*num_bowties, 1);
for i =1:4*num_bowties
    if mod(i,4) ==1
        if mod(i,4*num_in_col) ==1
            f(i) = 1/Cwire *(-x(i)/Rwire - (x(i) -x(i+1))/Re -(x(i)-x(i+4))/Rwire + u.Wire1Bias);
        elseif mod(i,4*num_in_col) == (4*num_in_col-3)
            f(i) = 1/Cwire *((x(i-4)-x(i))/Rwire - (x(i) -x(i+1))/Re);
        else
            f(i) = 1/Cwire *((x(i-4)-x(i))/Rwire - (x(i) -x(i+1))/Re -(x(i)-x(i+4))/Rwire);
        end
    elseif mod(i,4) ==2
    %f(i) = ((-x(i)/Re - A*Jnano(phi,beta*(x(i) - x(i+1))/d, ROC, taby) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i+1)/Rc +A*Jnano(phi,beta*(x(i) - x(i+1))/d, ROC, taby)+u(i+1)))/(1-Cec^2/(Cec+Cp)^2);
    %disp(Jnano(phi,beta*(x(i) - x(i+1))/d, ROC, taby));
    %disp(Jnano(phi,beta*(x(i) - x(i+1))/d, ROC, taby));
    f(i) = (((x(i-1)-x(i))/Re - Jn )/(Cec+Cp) +(Cec/(Cec+Cp)^2)*((x(i+2)-x(i+1))/Rc + Jn))/(1-Cec^2/(Cec+Cp)^2);
    elseif mod(i,4) ==3
    %f(i) = ((-x(i)/Rc + A*Jnano(phi,beta*(x(i-1) - x(i))/d, ROC, taby) +u(i))/(Cec+Cp) +(Cec/(Cec+Cp)^2)*(-x(i-1)/Re -A*Jnano(phi,beta*(x(i-1) - x(i))/d, ROC, taby)+u(i-1)))/(1-Cec^2/(Cec+Cp)^2);   
    f(i) = (((x(i+1)-x(i))/Rc + Jn)/(Cec+Cp) +(Cec/(Cec+Cp)^2)*((x(i-2)-x(i-1))/Re - Jn ))/(1-Cec^2/(Cec+Cp)^2);
    else
        if mod(i, 4*num_in_col) ==4
            f(i) =1/Cwire*(-x(i)/Rwire + (x(i) -x(i-1))/Rc -(x(i)-x(i+4))/Rwire + u.Wire2Bias);
        elseif mod(i, 4*num_in_col) ==0
            f(i) =1/Cwire*((x(i-4)-x(i))/Rwire + (x(i) -x(i-1))/Rc);
        else
        f(i) =1/Cwire*((x(i-4)-x(i))/Rwire + (x(i) -x(i-1))/Rc -(x(i)-x(i+4))/Rwire);
        end 
    end 
end 