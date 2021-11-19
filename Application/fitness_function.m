function [fitness_value] = fitness_function(centroids)

    fishertable = readtable('fisheriris.csv');
    X = table2array(fishertable(1:10,1:2));
    
   [num_samples,~] = size(X);
   minimums = zeros(num_samples,1);
   for i=1:10
        temp = zeros(1,2);
        c=1;
        for j=1:4
            if rem(j,2)==0
               j=j+1;
               continue
            end
        temp(c)=sqrt(sum(power((X(i,:)-centroids(j:j+1)),2)));
        c=c+1;
        end
            minimums(i) = min(temp);    
    end
    fitness_value = sum(minimums);
    
end


