function f = f_special(x0, p, u)
A = zeros(length(x0));
row = sqrt(length(x0));
col = sqrt(length(x0));
for i = 1:length(x0)
    A(i,i) = 7;
    if mod(i,2) == 1
        A(i,i+1)= A(i,i+1)-5;
        if  i> row
            A(i,i) = A(i,i) + 0.5;
            A(i,i-row+1) = A(i,i-row+1) -0.5;
        end 
    else
        A(i,i-1)= -5;
        if i<= row*(col-1)
            A(i,i) = A(i,i) + 0.5;
            A(i,i+row-1) = A(i,i+row-1) - 0.5;
        end
    end
end
f = A*x0;
end