function power = covariate_net_reg(subject_num, N, sigma0,...
    f2, DF, rho_in, rho_out, cluster_size, Covariance, ...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel, Method)
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

%   Check the number of input arguments and set default method if needed
if nargin == 14  % If there are 14 input arguments
    Method = 'BNPower';  % Set the default method to 'BNPower'
end

% Setup the Total Iterations and Waitbar:
wb = waitbar(0, 'Please wait...');

% Disable warnings
warning('off','all');

% Check for Parallel Computing Toolbox:
hasParallelToolbox = license('test','Distrib_Computing_Toolbox');

if hasParallelToolbox && Useparallel
    % Check for Parallel Pool:
    if isempty(gcp('nocreate'))
        % Get the default cluster profile
        cluster = parcluster('local');

        % Get the maximum number of available workers
        numWorkers = cluster.NumWorkers;

        % Start the parallel pool with all available workers
        parpool('local', numWorkers);
    end

    % Initialize Data Queue:
    dq = parallel.pool.DataQueue;
    % Setup Listener for Progress Updates:
    afterEach(dq, @(varargin) iIncrementWaitbar(wb));
else
    % No data queue for serial execution
    dq = [];
end

power_list = zeros(M_rep,1);

% Main power calculation loop
if hasParallelToolbox && Useparallel
    L = chol(Covariance);
    parfor m = 1:M_rep

        % Initialize Clist and CID at the start of each iteration
        Clist = [];  % or some appropriate initial value
        CID = [];    % or some appropriate initial value
        thresh = 0.8;
        
        % Perform sampling and network detection using greedy algorithm
        [W1,Wp,Clist_GT,Wm, X, threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, DF, cluster_size, rho_in, rho_out, L, Reliability);
        
        % select method
        if strcmp(Method,'BNPower')
            % greedy
            func = @(C) sum(sum(C))/size(C,1);
            [~, Clist, CID] = greedy(Wp, func, threshold_GT);
        elseif strcmp(Method,'NBS')
            thresh = 0.5;
            % NBS for network analysis
            [~, Clist, CID] = NBS(Wp, -log10(threshold_GT));
        end
        
        % Filter based on Jaccard index to proceed with permutation test
        J = length(intersect(Clist_GT(1:cluster_size),Clist(1:CID)))/length(union(Clist_GT(1:cluster_size),Clist(1:CID)));
        if J >= thresh
            P_value = permutation_test_reg(Wm, X, Wp, CID, Clist, M_perm, threshold_GT, Method);
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
        thresh = 0.8;
        [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, DF, cluster_size, rho_in, rho_out, L, Reliability);
        % select method
        if strcmp(Method,'BNPower')
            % greedy
            func = @(C) sum(sum(C))/size(C,1);
            [~, Clist, CID] = greedy(Wp, func, threshold_GT);
        elseif strcmp(Method,'NBS')
            thresh = 0.5;
            % NBS for network analysis
            [~, Clist, CID] = NBS(Wp, -log10(threshold_GT));
        end

        % filter by Jaccard index
        J = length(intersect(Clist_GT(1:cluster_size),Clist(1:CID)))/length(union(Clist_GT(1:cluster_size),Clist(1:CID)));
        if J >= thresh
            P_value = permutation_test_reg(Wm, X, Wp, CID, Clist, M_perm, threshold_GT, Method);
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