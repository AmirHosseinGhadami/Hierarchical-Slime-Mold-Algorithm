function [InitFunction, CostFunction, FeasibleFunction] = Sphere

InitFunction = @SphereInit;
CostFunction = @SphereCost;
FeasibleFunction = @SphereFeasible;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = SphereInit(OPTIONS,trial)

global MinParValue MaxParValue
load Pop_Sphere_100_30 pop  
Granularity = 0.1;
MinParValue = -100*ones(1,OPTIONS.numVar);
MaxParValue = 100*ones(1,OPTIONS.numVar);
% Initialize population
if trial==-1
    for popindex = 1 : 1
    chrom = MinParValue + (MaxParValue - MinParValue + 1) .* rand(1,OPTIONS.numVar);
    Population.chrom = chrom;
end
else
for popindex = 1 : OPTIONS.popsize
     chrom=pop(trial).Population(popindex,:);
     Population(popindex).chrom = chrom;
%     Population(popindex).pBestChrom = chrom;

end
end

% Cost= sphere1Cost(OPTIONS,Population)
% for j=  1:OPTIONS.popsize  
% Population(j).pBestCost=Cost(j)
% end

OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SphereCost(OPTIONS, Population)

% Compute the cost of each member in Population
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    sum = 0;   
    for i = 1 : p
        x = Population(popindex).chrom(i);
        sum = sum + x^2;
    end
    Population(popindex).cost = sum;    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SphereFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue(k));
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue(k));
    end
end
return;