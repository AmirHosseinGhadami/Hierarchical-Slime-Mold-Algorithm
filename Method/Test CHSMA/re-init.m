function h =init(X,AllFitness, N, dim, h_number, it, lb, ub)

[SmellOrder,SmellIndex] = sort(AllFitness); 

lu = ones(2,dim);
lu(1,:) = lu(1,:)*(lb(1,1));
lu(2,:) = lu(2,:)*(ub(1,1));
h.rank = ones(1,h_number);
for n = 1:h_number  
    
    h_len = (N-h_number)/h_number;
    start_range = h_number+ h_len*(n-1);
    m = zeros(1,h_len);
    m(1,1)= n;
    for counter = 2:h_len
        m (1, counter) = start_range; 
        start_range = start_range +1 ;
    end
    
    for i = 1:h_len
        h(n).X(i,:) = X(SmellIndex(m(1,i)),:);
        h(n).F(i,1) = SmellOrder(m(1,i),1);
        h(n).pBest_X(i,:) = X(SmellIndex(i),:);
        h(n).pBest_F(i,1) = SmellOrder(i,1);
        h(n).pre_F(i,1) = SmellOrder(m(1,i),1);
        
        % initialize smell path archive
        r_temp = rand(1,dim);
        m_temp = r_temp>0.5; 
        N_temp = -1 * ones(1,dim);
        h(n).Levy_Archive(i,:) = (N_temp+m_temp)+ m_temp;

        h(n).Levy_Change_Order(i,:) = randperm(dim);
        h(n).Levy_Change_Counter(i,1) = 1; 
        h(n).path_flag(i,1) = 1 ;

    end
        %% Defining each hierarchy's CMA's parameters  

%         CMA.xmeanw = (lu(1, :) + rand(1, dim) .* (lu(2, :) - lu(1, :)))'; % object parameter start point
        CMA.xmeanw =  mean(h(n).X(:,:))';
%         CMA.xmeanw =  mean(h(n).X(1,:))';
        CMA.sigma = 0.25; %0.25; % CMA.sigma = 0.25;
        CMA.minsigma = 1e-5;   %15;
        CMA.maxsigma = max(lu(2, :)-lu(1, :)) / sqrt(dim); % initial step size, minimal step size
        CMA.flginiphase = 1; % Initial phase

        % Parameter setting: selection,
%         CMA.lambda = 4 + floor(3 * log(dim)); CMA.mu = floor(CMA.lambda/2);
        CMA.lambda = 20; %CMA.mu = floor(CMA.lambda/2);
        CMA.mu = 10 ;  % 10; for 5 pop

%         CMA.lambda = 4 + floor(2 * log(dim)); CMA.mu = floor(CMA.lambda/3);
%         CMA.lambda = 30; CMA.mu = floor(CMA.lambda/7);

        CMA.arweights = log((CMA.lambda + 1)/2) - log(1:CMA.mu)'; % muXone array for weighted recomb.
        % lambda = 10; mu = 2; arweights = ones(mu, 1); % uncomment for (2_I, 10)-ES
        % parameter setting: adaptation
        CMA.cc = 4/(dim + 4); CMA.ccov = 2/(dim + 2^0.5)^2;
        CMA.cs = 4/(dim + 4); CMA.damp = (1 - min(0.7, dim * CMA.lambda/it)) / CMA.cs + 1;

        % Initialize dynamic strategy parameters and constants
        CMA.B = eye(dim); CMA.D = eye(dim); CMA.BD = CMA.B * CMA.D; CMA.C = CMA.BD * transpose(CMA.BD);
        CMA.pc = zeros(dim, 1); CMA.ps = zeros(dim, 1);
        CMA.cw = sum(CMA.arweights) / norm(CMA.arweights); CMA.chiN = dim^0.5 * (1 - 1/(4 * dim) + 1/(21 * dim^2));

        % Generation loop
        %     disp(['  (' num2str(mu) ', ' num2str(lambda) ')-CMA-ES (w = [' num2str(arweights', '%5.2f') '])' ]);
%         counteval = 0; flgstop = 0;
%         mincost=[];
        % Boundary
        CMA.lb = (ones(CMA.lambda, 1) * lu(1, :))';
        CMA.ub = (ones(CMA.lambda, 1) * lu(2, :))';
        
        
        
        
   
    h(n).CMA = CMA;
    
end
end
