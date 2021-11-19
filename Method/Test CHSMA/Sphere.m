function [F] = Sphere(X)

popsize = size(X,1);
F= zeros(popsize);
p = size(X,2);
for popindex = 1 : popsize
    sum = 0;   
    for i = 1 : p
        x = X(i);
        sum = sum + x^2;
    end
    F(popindex) = sum;    
end
return

