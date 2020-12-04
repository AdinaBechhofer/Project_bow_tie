function f = f_linear(x, t,p,u,b)
C = zeros(length(x));
G = zeros(length(x));
% Stamp C and invert; stamp G
row = sqrt(length(x));
col = row;
for i = 1:length(x)
    C(i,i) = 5+ 0.5;
    if mod(i,2) == 1
        C(i,i+1)= C(i,i+1)-5;
        G(i,i) = G(i,i) +1/5.8;
        if  i> row
            C(i,i) = C(i,i) + 0.7;
            C(i,i-row+1) = C(i,i-row+1) -0.7;
        end 
    else
        C(i,i-1)= C(i,i-1) -5;
        G(i,i) = 1/5.8;
        if i<= row*(col-1)
            C(i,i) = C(i,i) +0.7;
            C(i,i+row-1) = C(i,i+row-1) - 0.7;
        end
    end
end

B = zeros(length(x), 2);
B(1:2:end, 1) = 1;
B(2:2:end, 2) =1;
f = G*C*x + B*u;
end 