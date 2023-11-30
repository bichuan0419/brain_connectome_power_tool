function power = covariate_net_reg(subject_num, N, sigma0,...
    f2, DF, rho_in, rho_out, cluster_size, Covariance, ...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel)
%COVARIATE_NET_REG Estimates power of network based on covariate regression.
%
%   This function estimates the statistical power of detecting a network
%   change when using a covariate regression approach on whole-brain 
%   connectome data. It uses a permutation testing framework to assess 
%   significance and applies a greedy algorithm for network identification.
%
%   Inputs:
%       subject_num   - Number of subjects in the dataset.
%       N             - Number of nodes in the network.
%       sigma0        - Baseline standard deviation for FC.
%       f2            - Cohen's f2 of the edges within the covariate-related 
%                       subnetwork.
%       DF            - Degrees of freedom.
%       rho_in        - Proportion of significant edges within the 
%                       covariate-related subnetwork.
%       rho_out       - Proportion of significant edges outside the
%                       covariate-related subnetwork.
%       cluster_size  - Size of the clusters in the network.
%       Covariance    - Covariance matrix for the network nodes.
%       FWER_threshold- Family-wise error rate threshold for significance.
%       M_rep         - Number of repetitions (or number of datasets)
%                       for power calculation.
%       M_perm        - Number of permutations for the permutation test.
%       Reliability   - Measure of test-retest reliability.
%
%   Output:
%       power         - Estimated statistical power for the network test.
%
%   See also PERMUTATION_TEST_REG, SAMPLING_IND_REG.


% Setup the Total Iterations and Waitbar:
wb = waitbar(0, 'Please wait...');

% Check for Parallel Computing Toolbox:
hasParallelToolbox = license('test','Distrib_Computing_Toolbox');

% Parallel pool setup or serial execution based on toolbox availability
if hasParallelToolbox && Useparallel
    % Initialize the parallel pool if not already started
    if isempty(gcp('nocreate'))
        cluster = parcluster('local');
        numWorkers = cluster.NumWorkers;
        parpool('local', numWorkers);
    end
    % Initialize Data Queue for progress updates in parallel execution
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(varargin) iIncrementWaitbar(wb));
else
    % Data Queue is not required for serial execution
    dq = [];
end

% Initialize list to hold power calculation results for each repetition
power_list = zeros(M_rep,1);

% Main power calculation loop
if hasParallelToolbox && Useparallel
    L = chol(Covariance);
    parfor m = 1:M_rep
        % Perform sampling and network detection using greedy algorithm
        [W1,Wp,Clist_GT,Wm, X, threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, DF, cluster_size, rho_in, rho_out, L, Reliability);
        func = @(C) sum(sum(C))/length(C); % Define function handle for greedy algorithm
        [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func, threshold_GT);
        
        % Filter based on Jaccard index to proceed with permutation test
        J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
        if J >= 0.8
            P_value = permutation_test_reg(Wm, X, Wp, CID_greedy, Clist_greedy, M_perm, func, threshold_GT);
            power_list(m) = double(P_value <= FWER_threshold);
        else
            power_list(m) = 0;
        end
        send(dq, m); % Update progress for parallel execution
    end
else
    % Serial execution if Parallel Computing Toolbox is not available
    L = chol(Covariance);
    for m = 1:M_rep
        [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, DF, cluster_size, rho_in, rho_out, L, Reliability);
        % greedy
        func = @(C) sum(sum(C))/length(C);
        
        [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func,threshold_GT);

        % filter by Jaccard index
        J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
        if J >= 0.8
            P_value = permutation_test_reg(Wm, X, Wp, CID_greedy, Clist_greedy, M_perm, func,threshold_GT);
            power_list(m) = double(P_value<= FWER_threshold);
            %         disp(double(P_value<= FWER_threshold))
        else
            power_list(m) = 0;
        end
        % Directly update the waitbar for serial execution
        waitbar(m/M_rep, wb, sprintf('Please wait...'));


    end
end
power = mean(power_list);
% Close the Waitbar:
close(wb);
end