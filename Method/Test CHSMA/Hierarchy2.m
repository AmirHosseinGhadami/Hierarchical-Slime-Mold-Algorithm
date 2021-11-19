% Source codes demo version 1.0
% ------------------------------------------------------------------------------------------------------------
% Main paper (Please refer to the main paper):
% Slime Mould Algorithm: A New Method for Stochastic Optimization
% Shimin Li, Huiling Chen, Mingjing Wang, Ali Asghar Heidari, Seyedali Mirjalili
% Future Generation Computer Systems,2020
% DOI: https://doi.org/10.1016/j.future.2020.03.055
% https://www.sciencedirect.com/science/article/pii/S0167739X19320941
% ------------------------------------------------------------------------------------------------------------
% Website of SMA: http://www.alimirjalili.com/SMA.html
% You can find and run the SMA code online at http://www.alimirjalili.com/SMA.html

% You can find the SMA paper at https://doi.org/10.1016/j.future.2020.03.055
% Please follow the paper for related updates in researchgate: https://www.researchgate.net/publication/340431861_Slime_mould_algorithm_A_new_method_for_stochastic_optimization
% ------------------------------------------------------------------------------------------------------------
%  Main idea: Shimin Li
%  Author and programmer: Shimin Li,Ali Asghar Heidari,Huiling Chen
%  e-Mail: simonlishimin@foxmail.com
% ------------------------------------------------------------------------------------------------------------
%  Co-author:
%             Huiling Chen(chenhuiling.jlu@gmail.com)
%             Mingjing Wang(wangmingjing.style@gmail.com)
%             Ali Asghar Heidari(aliasghar68@gmail.com, as_heidari@ut.ac.ir)
%             Seyedali Mirjalili(ali.mirjalili@gmail.com)
%             
%             Researchgate: Ali Asghar Heidari https://www.researchgate.net/profile/Ali_Asghar_Heidari
%             Researchgate: Seyedali Mirjalili https://www.researchgate.net/profile/Seyedali_Mirjalili
%             Researchgate: Huiling Chen https://www.researchgate.net/profile/Huiling_Chen
% ------------------------------------------------------------------------------------------------------------
% _____________________________________________________
%  Co-author and Advisor: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%
%       Homepage: http://www.alimirjalili.com
% _____________________________________________________
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Max_iter: maximum iterations, N: populatoin size, Convergence_curve: Convergence curve
% To run SMA: [Destination_fitness,bestPositions,Convergence_curve]=SMA(N,Max_iter,lb,ub,dim,fobj)
function [Destination_fitness,bestPositions,Convergence_curve]=Hierarchy2(N,Max_iter,lb,ub,dim,fobj)
% disp('SMA is now tackling your problem')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_name = 'cec17_func';
fitness_func = str2func(func_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func_name = 'Sphere';
% fitness_func = str2func(func_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
beta = 1.4;
%%
Number_of_runs = 1 ;
Convergence_curve=zeros(Number_of_runs,Max_iter);  % Change this if you need less runs
%%
for Function_number=1:30
    disp(append('Function ',string(Function_number)));
    
Function_name = append('F',string(Function_number));
    
for count_run = 1:Number_of_runs % to run 30 individual times

% initialize position
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb);


% X = load(append('Populations\Pop_',Function_name,'_P100_D',string(dim))).Positions; 
% X = X(1:dim,:);
% X = load(append('input_data\M_',Function_name,'_D',string(dim))).Positions; 
% X = load(append('input_data_',func_name,'\M_',string(Function_number),'_D100.txt')); 
% X = X(1:N,1:dim); %  M_1_D10
X = load(append('Population/p100_d100_f',string(Function_number))); 
X = X.Positions(1:N,1:dim); %  M_1_D10

it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%         AllFitness(i) = fobj(X(i,:));
        AllFitness(i)=feval(fitness_func,X(i,:)',Function_number);
%         AllFitness(i) = feval(fitness_func,X(i,:));
    end
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    S = bestFitness-worstFitness+eps;
        Convergence_curve(count_run,1)=bestFitness;
        Destination_fitness = bestFitness;
        
        h_number = 5; 
        h_len = (N)/h_number;
        h =init(X,AllFitness, N,dim,h_number,it,lb,ub);
%%


% Main loop
while  it <= Max_iter -1 
    
    %%
    tempX = [] ;
    tempF = [];
    tempBest = [] ; 
    for n = 1:h_number
        tempX = [tempX; h(n).X];
        tempF = [tempF; h(n).F];
        tempBest = [tempBest; h(n).pBest_X];
    end
    
    [SmellOrder,SmellIndex] = sort(tempF);  %Eq.(2.6)
    temp = zeros(N,dim);
    tempBest = zeros(N,dim);
    for i = 1:N
        temp(i,:) = tempX(SmellIndex(i),:);
        tempBest(i,:) = tempBest(SmellIndex(i),:);
    end
  %%
    
    
    for n = 1:h_number
                t_h = h;
                h(n).CMA.arz = randn(dim, h(n).CMA.lambda);
    %             h(n).CMA.arx = h(n).CMA.xmeanw * ones(1, h(n).CMA.lambda) + h(n).CMA.sigma * (h(n).CMA.BD * h(n).CMA.arz);
                h(n).CMA.arx = h(n).CMA.xmeanw * ones(1, h(n).CMA.lambda) + h(n).CMA.sigma * (h(n).CMA.BD * h(n).CMA.arz);
%                 h(n).CMA.arx = h(n).X(1,:)' * ones(1, h(n).CMA.lambda) + h(n).CMA.sigma * (h(n).CMA.BD * h(n).CMA.arz);

                               % h(n).X(1,:)

                % Handle the elements of the variable which violate the boundary
                I = find(h(n).CMA.arx > h(n).CMA.ub);
                h(n).CMA.arx(I) = 2 * h(n).CMA.ub(I) - h(n).CMA.arx(I);
                aa = find(h(n).CMA.arx(I) < h(n).CMA.lb(I));
                h(n).CMA.arx(I(aa)) = h(n).CMA.lb(I(aa));
                I = find(h(n).CMA.arx < h(n).CMA.lb);
                h(n).CMA.arx(I) = 2 * h(n).CMA.lb(I) - h(n).CMA.arx(I);
                aa = find(h(n).CMA.arx(I) > h(n).CMA.ub(I));
                h(n).CMA.arx(I(aa)) = h(n).CMA.ub(I(aa));
                
                U = h(n).CMA.arx';
                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 h(n).CMA.arfitness = feval(fitness_func,U',Function_number);%benchmark_func(U, problem, o, A, M, a, alpha, b);
%                 h(n).CMA.arfitness = feval(fitness_func,U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%                 h(n).CMA.arfitness=h(n).CMA.arfitness';
%                 counteval = counteval + lambda;
%                 [h(n).CMA.arfitness, h(n).CMA.arindex] = sort(h(n).CMA.arfitness); % minimization
                h(n).CMA.xold = h(n).CMA.xmeanw; % for speed up of Eq. (14)

%                 h(n).X(1:h_len,:) = h(n).CMA.arx(:,h(n).CMA.arindex(1:h_len))';


                
                if ((h(1).rank(1,n)) == 1)
                    t_position = bestPositions;
                else
                    rank_indexes = (find(h(1).rank<h(1).rank(1,n)));
                    t_position = h(rank_indexes(randperm(length(rank_indexes),1))).pBest_X(randperm(1,1),:); 
                end
                
                
                a = tanh(it/Max_iter);
                b = tanh(1-it/Max_iter);
                c = (1-it/Max_iter);
                
                S = abs(bestFitness-worstFitness+eps);
                factor =  abs(bestFitness(1)-h(n).F(1))/S;
                
                alpha = it/Max_iter;
                leader = (factor) * h(n).X(1,:) + (1-factor) * h(n).CMA.arx(:,1)';
%                 leader = (factor) * h(n).X(1,:) + (1-factor) * h(n).CMA.arx(:,1)';
                %%
                if  it>2
                    union = [h(n).X ; h(n).archive.X];
                    union = union(randperm(size(union,1),1), :);
                else
                    union = bestPositions;
                end
                %%
                if rand()< 1  %> Max_iter/2
                    h(n).X(1,:) = h(n).CMA.arx(:,1)' +   rand()*(bestPositions - h(n).X(1,:));% + b*rand()*(-t_position - h(n).X(1,:)) ; % + c* rand()* (h(n).X(r1,:)-union);
%                     h(n).X(1,:) = leader +   factor*rand()*(t_position - h(n).X(1,:)); % + h(n).Factor(1,1)*(h(n).X(r1,:)-union) ; %   +  b*sum(levy(5,dim,beta)) ; %+ (it/Max_iter)*rand()*(bestPositions(1,:)- h(n).X(1,:));
                else                 %  rand()*b*(t_position - h(n).X(1,:))+
%                     r1 = randperm(h_len,1);
                    h(n).X(1,:) = h(n).CMA.arx(:,1)' +   h(n).CMA.sigma*rand()*(t_position - h(n).X(1,:));% + c* rand()*(tempBest(randperm(25,1),:)- h(n).X(1,:)) ; % + h(n).Factor(1,1)*(h(n).X(r1,:)-union) ; %   +  b*sum(levy(5,dim,beta)) ; %+ (it/Max_iter)*rand()*(bestPositions(1,:)- h(n).X(1,:));
%                     h(n).X(1,:) = leader + rand()*(tempBest(randperm(25,1),:)- h(n).X(1,:));
                end                                  


for i = 2: h_len
%     S = abs(bestFitness-worstFitness+eps);
%     factor =  abs(bestFitness(1)-h(n).F(i))/S
    
    
    if  it>2
%         randperm(size(h(n).archive,1),1);
        union = [h(n).X ; h(n).archive.X];
        union = union(randperm(size(union,1),1), :);
    else
        union = h(n).X(i,:);
    end
                 
    S = h(n).F(1)-h(n).F(h_len)+eps;
    factor = abs((h(n).F(1) - h(n).F(2)))/S;

    follower = (factor) * h(n).X(i,:) + (1-factor) * h(n).CMA.arx(:,i)';

    
    if rand() < factor %factor * (1-it/Max_iter)  %|| h(1).rank(1,n) ==1 
        r1 = randperm(round(h_len/2),1);
%         h(n).X(i,:) = follower  +  c* rand()*(tempX(randperm(25,1),:)- h(n).X(i,:));% + h(n).Factor(i,1)*(h(n).X(r1,:)-union) ;% + rand()*(h(n).X(r1,:)-union);% + c*rand()*(h(h_n).pBest_X(r_i(1,1),:)-union);% + b* rand()*(h(h_n).pBest_X(r_i(1,2),:)-union) ;  % +  a* rand()*(h(n).X(1,:)-t_h(n).X(i,:)) 
        h(n).X(i,:) = follower  +  c* factor *(h(n).pBest_X(randperm(i-1,1),:)- h(n).X(i,:));% + h(n).Factor(i,1)*(h(n).X(r1,:)-union) ;% + rand()*(h(n).X(r1,:)-union);% + c*rand()*(h(h_n).pBest_X(r_i(1,1),:)-union);% + b* rand()*(h(h_n).pBest_X(r_i(1,2),:)-union) ;  % +  a* rand()*(h(n).X(1,:)-t_h(n).X(i,:)) 
    else

        r1 = randperm(h_len,1);
        h(n).X(i,:) = h(n).CMA.arx(:,i)';    
    end
end
  



    end
    
    
    
    
    [h,bestPositions,Destination_fitness,bestFitness,worstFitness] = h_sort (h, h_number, bestPositions, Destination_fitness, ub , lb, dim,fitness_func,Function_number);

    for n = 1:h_number
         h(n) = adapt(h(n),dim,it, Max_iter ,bestPositions,Destination_fitness);
    end
    
    
    
    
    Convergence_curve(count_run,it+1)=Destination_fitness;
    disp(string(Destination_fitness));
    it=it+1;
end

end % end count run loop

    disp(string(Destination_fitness));

% name = append('Results\SMA_',fobj,'_',Function_name);
% save(name, 'Convergence_curve');% Number

% display(['Average Cost is: ',num2str(mean(CostMatrix))])
end
end


function [z] = levy(n,m,beta)
% This function implements Levy's flight. 
% For more information see 
%'Multiobjective cuckoo search for design optimization Xin-She Yang, Suash Deb'. 
% Coded by Hemanth Manjunatha on Nov 13 2015.
% Input parameters
% n     -> Number of steps 
% m     -> Number of Dimensions 
% beta  -> Power law index  % Note: 1 < beta < 2
% Output 
% z     -> 'n' levy steps in 'm' dimension
    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator
    sigma_u = (num/den)^(1/beta);% Standard deviation
    u = random('Normal',0,sigma_u^2,n,m); 
    
    v = random('Normal',0,1,n,m);
%     z = u./(abs(v).^(1/beta));  % The orignial one
    
    z = abs(u)./(abs(v).^(1/beta));    % We use this one to apply the archive
    % Just the "u" is in abs funciton to only get the positive value
     z = (u)./((v).^(1/beta));  % abs
     
end
