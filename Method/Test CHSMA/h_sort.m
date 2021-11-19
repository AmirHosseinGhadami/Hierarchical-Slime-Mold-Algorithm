function [h,bestPositions,Destination_fitness,bestFitness,worstFitness] = h_sort (h, h_number, bestPositions, Destination_fitness, ub , lb, dim,fitness_func,Function_number)
bestFitness = 10^100;
worstFitness = -10^100;
for n = 1: h_number  
    l = size(h(n).X,1);    
    t=h(n);
    X= t.X;
    counter = 1;
    
    archive_candidates = [];
    
    for i=1:l
        % Check if solutions go outside the search space and bring them back
        Flag4ub=t.X(i,:)>ub;
        Flag4lb=t.X(i,:)<lb;
        t_genes = (ub-lb).*rand(1,dim) + lb;
        t.X(i,:)=(t.X(i,:).*(~(Flag4ub+Flag4lb)))+t_genes.*Flag4ub+t_genes.*Flag4lb;
%         AllFitness(i) = fobj(X(i,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AllFitness(i)=feval(fitness_func,X(i,:)',Function_number);
%         AllFitness(i)=feval(fitness_func,X(i,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    end
    
    %%

       
    
%%
    [valParents, I] = min([h(n).pre_F(:), AllFitness'], [], 2);
%     t.X(I==1, :) = t.pre_X(I==1,:);   
    
    h(n).archive = updateArchive(h(n).archive, t.pre_X(I == 2, :), h(n).pre_F(I == 2));
    
    h(n).goodCR = h(n).Factor2(I == 2);
    h(n).goodF = h(n).Factor(I == 2);
    [h(n).Factor , h(n).Factor2] = randFCR(l, h(n).Fm2, 0.1, h(n).Fm, 0.1);
%%    
    

    [SmellOrder,SmellIndex] = sort(valParents);  %Eq.(2.6)
    
%     bestFitness = SmellOrder(1);
    if bestFitness > SmellOrder(1)
        bestFitness = SmellOrder(1);
    end
    if  SmellOrder(h_number) > worstFitness
        worstFitness = SmellOrder(h_number);
    end

    for i=1:l           
        
        h(n).X(i,:) = t.X(SmellIndex(i),:);
        h(n).F(i,1) = SmellOrder(i);
        %% Updating the archive candidates
%         if i< 4  && (t.pre_F(SmellIndex(i))-h(n).F(i,1))>0
%             archive_candidates(counter,:) = t.pre_X(SmellIndex(i),:);
%             counter = counter+1;
%         end
        %%
        h(n).pre_F(i,1) = t.F(SmellIndex(i),1);
        h(n).pre_X(i,:) = t.X(SmellIndex(i),1);
        
        if SmellOrder(i)< t.pBest_F(SmellIndex(i),:)            
            h(n).pBest_X(i,:) = t.X(SmellIndex(i),:);
            h(n).pBest_F(i,1) = SmellOrder(i);
        else
            h(n).pBest_X(i,:) = t.pBest_X(SmellIndex(i),:);
            h(n).pBest_F(i,1) = t.pBest_F(SmellIndex(i),1);
        end        
        h(n).Levy_Archive(i,:) = t.Levy_Archive(SmellIndex(i),:);
        h(n).Levy_Change_Order(i,:) = t.Levy_Change_Order(SmellIndex(i),:);
        h(n).Levy_Change_Counter(i,1) = t.Levy_Change_Counter(SmellIndex(i),1); 
        h(n).path_flag(i,1) = t.path_flag(SmellIndex(i),1);                
    end   
    

        %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   if size(archive_candidates,1) ~= 0   
%       
%       index_number = randi(ceil(size(archive_candidates,1)/2));
%       indexes = randperm(size(archive_candidates,1), index_number);
%         
%       for ii=1:index_number
%           ar_len = round(dim/3);
%           if  size(h(n).archive,1) >= ar_len
%               h(n).archive(randperm(ar_len,1),:) = archive_candidates(indexes(ii),:);
%           else
%               arch_len = length(h(n).archive);
%               h(n).archive(arch_len+ii,:) = archive_candidates(indexes(ii),:);
%           end
%       end
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
  
%%% For ranking Hierarchies
  leaders= [];
  leaders_Fitness =[];
  for n=1:h_number
      leaders= [leaders; h(n).X(1,:)];
      leaders_Fitness = [leaders_Fitness; h(n).F(1,1)];
  end
  [n_,leaderIndex] = sort(leaders_Fitness); 
  h(1).rank(1,:) = leaderIndex; 
 
  

        
          
          
%   if length(h(n).archive)< l-index
%       del_indexes = randperm(l,index_number);
%       for m=1:index_number
%           h(n).archive(del_indexes(m),:) = archive_candidates(m,:);
%       end
%   else
%       h(n).archive(del
%       
%   end
end