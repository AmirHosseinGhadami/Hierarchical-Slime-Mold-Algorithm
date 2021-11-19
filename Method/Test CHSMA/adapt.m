function h = adapt(h,dim,it, Max_iter, bestPositions,Destination_fitness)
    
    CMA= h.CMA;
    ub = CMA.ub;
    lb = CMA.lb;
    
    
    

%     counteval = counteval + lambda;
    % Sort by fitness and compute weighted mean in xmeanw
%     CMA.arfitness = h.F; % minimization
%     CMA.xold = CMA.xmeanw; % for speed up of Eq. (14)
    
%     temp_x = [bestPositions(1,:);h.X(1:CMA.mu-1,:)];

    
    CMA.arfitness=h.F';
%                 counteval = counteval + lambda;
    [CMA.arfitness, CMA.arindex] = sort(CMA.arfitness); % minimization
    
    temp_x = h.X(1:CMA.mu,:);
    CMA.xmeanw = temp_x' * CMA.arweights/sum(CMA.arweights);
    CMA.zmeanw = CMA.arz(:, CMA.arindex(1:CMA.mu)) * CMA.arweights/sum(CMA.arweights);
    
%     CMA.xmeanw = h.X(1:CMA.mu,:)' * CMA.arweights/sum(CMA.arweights);
%     CMA.zmeanw = h.X(1:CMA.mu,:)' * CMA.arweights/sum(CMA.arweights);
    
%     h(n).CMA.xmeanw = h(n).CMA.arx(:, h(n).CMA.arindex(1:h(n).CMA.mu)) * h(n).CMA.arweights/sum(h(n).CMA.arweights);
%     h(n).CMA.zmeanw = h(n).CMA.arz(:, h(n).CMA.arindex(1:h(n).CMA.mu)) * h(n).CMA.arweights/sum(h(n).CMA.arweights);
    
    % Adapt covariance matrix
    CMA.pc = (1-CMA.cc) * CMA.pc + (sqrt(CMA.cc * (2-CMA.cc)) * CMA.cw/CMA.sigma) * (CMA.xmeanw-CMA.xold); % Eq. (14)
    if ~CMA.flginiphase % do not adapt in the initial phase
        CMA.C = (1-CMA.ccov) * CMA.C + CMA.ccov * CMA.pc * transpose(CMA.pc);           % Eq. (15)
    end
    % adapt sigma
    CMA.ps = (1-CMA.cs) * CMA.ps + (sqrt(CMA.cs * (2-CMA.cs)) * CMA.cw) * (CMA.B * CMA.zmeanw);      % Eq. (16)
    
    CMA.damp = (1 - min(0.7, dim * CMA.lambda/it)) / CMA.cs + 1;
    CMA.sigma = CMA.sigma * exp((norm(CMA.ps)-CMA.chiN)/CMA.chiN/CMA.damp);        % Eq. (17)
%     CMA.sigma = CMA.sigma * exp((norm(CMA.ps)-CMA.chiN)/CMA.chiN);        % Eq. (17)


    % Update B and D from C
    if mod(it/CMA.lambda, 1/CMA.ccov/dim/5) < 1
        CMA.C = triu(CMA.C) + transpose(triu(CMA.C, 1)); % enforce symmetry
        [CMA.B, CMA.D] = eig(CMA.C);
        % limit condition of C to 1e14 + 1
        if max(diag(CMA.D)) > 1e14 * min(diag(CMA.D))
            tmp = max(diag(CMA.D))/1e14 - min(diag(CMA.D));
            CMA.C = CMA.C + tmp * eye(dim); CMA.D = CMA.D + tmp * eye(dim);
        end
        CMA.D = diag(sqrt(diag(CMA.D))); % D contains standard deviations now
        CMA.BD = CMA.B * CMA.D; % for speed up only
    end % if mod

    % Adjust minimal step size
    if CMA.sigma * min(diag(CMA.D)) < CMA.minsigma ...
            | CMA.arfitness(1) == CMA.arfitness(min(CMA.mu + 1, CMA.lambda)) ...
            | CMA.xmeanw == CMA.xmeanw ...
            + 0.2 * CMA.sigma * CMA.BD(:, 1 + floor(mod(it/CMA.lambda, dim)))
        CMA.sigma = 1.4 * CMA.sigma;

        % flgstop = 1;
    end
    if CMA.sigma > CMA.maxsigma
        CMA.sigma = CMA.maxsigma;
    end

    % Test for end of initial phase
    if CMA.flginiphase & it/CMA.lambda > 2/CMA.cs
        if (norm(CMA.ps)-CMA.chiN)/CMA.chiN < 0.05 % step size is not much too small
            CMA.flginiphase = 0;
        end
    end
    h.CMA = CMA ; 
end