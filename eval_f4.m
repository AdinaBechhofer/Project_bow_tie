function f= eval_f4(x,t,p,u,b)
num_bowties = p.NumBowties;
%Re = p.REmitter;
%Rc = p.RCollector;
A = p.Area;
beta = p.Beta;
d = p.Distance;
phi = p.workFunction;
Cec = p.CemitterCollector;
Cp = p.Cparasitic;
ROC = p.Radius;
taby = p.taby;
Ivec = zeros(2*num_bowties,1);
C = zeros(2*p.NumBowties);
G = zeros(2*p.NumBowties); % this is for the linear part. Non-linear part is stamped in eval_f 

for i = 1:2*p.NumBowties
    C(i,i) = p.CemitterCollector+p.Cparasitic;
    if mod(i,2) == 1
        C(i,i+1)= C(i,i+1)-p.CemitterCollector;
        G(i,i) = G(i,i) -1./(p.REmitter);
        if  i> 2*p.row
            C(i,i) = C(i,i) + p.Ccoupling;
            C(i,i-2*p.row+1) = C(i,i-2*p.row+1) -p.Ccoupling;
        end 
    else
        C(i,i-1)= -p.CemitterCollector;
        G(i,i) = -1./(p.RCollector);
        if i<= 2*p.row*(p.col-1)
            G(i,i) = G(i,i) + p.Ccoupling;
            C(i,i+2*p.row-1) = C(i,i+2*p.row-1) + p.Ccoupling;
        end
    end
end

for i = 1:2*num_bowties
    if mod(i,2)==1
        Ivec(i) = - A*Jnano(phi,beta*(x(i)-x(i+1))/d,ROC,taby);
    else
        Ivec(i) = A*Jnano(phi,beta*(x(i-1)-x(i))/d,ROC,taby);
        %disp(A*Jnano(phi,beta*(x(i-1)-x(i))/d,ROC,taby))
    end
    
end
uvec = [u.vEmitter; u.vCollector];

f = C\(G*x+Ivec + b*uvec);
