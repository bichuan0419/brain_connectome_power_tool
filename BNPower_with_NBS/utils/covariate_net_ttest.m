function power = covariate_net_ttest(n1, n2, N, sigma0, ...
    Cohensd, rho_in, rho_out, cluster_size, Covariance,...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel, Method)

%COVARIATE_NET_TTEST Estimates power of data-driven network anlaysis based
% on a two-sample t-test.
%
%   This function estimates the statistical power of detecting a network
%   change using a two-sample t-test on whole-brain connectome data. It
%   uses a permutation testing framework to assess significance and applies
%   a greedy algorithm for network identification.
%
%   Inputs:
%       n1             - Number of subjects in group A.
%       n2             - Number of subjects in group B.
%       N              - Total number of nodes in the network.
%       sigma0         - Standard deviation of the edges.
%       Cohensd        - Effect size (Cohen's d).
%       rho_in         - Proportion of significant edges within the
%                        covariate-related subnetwork.
%       rho_out        - Proportion of significant edges outside the
%                        covariate-related subnetwork.
%       cluster_size   - Size of the clusters in the network.
%       Covariance     - Covariance matrix for the network nodes.
%       FWER_threshold - Family-wise error rate threshold for significance.
%       M_rep          - Number of repetitions for the simulation.
%       M_perm         - Number of permutations for the permutation test.
%       Reliability    - Measure of test-retest reliability.
%       Method         - Method to be used for the data-driven network-level
%                        power analysis.
%
%   Output:
%       power          - Estimated statistical power for the network test.

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

if hasParallelToolbox && Useparallel
    
    parfor m = 1:M_rep
         % Initialize Clist and CID at the start of each iteration
        Clist = [];  % or some appropriate initial value
        CID = [];    % or some appropriate initial value
        thresh = 0.8;

        theta = Cohensd*sigma0;
        if ~isempty(Covariance)
            L = chol(Covariance);
        else
            L = [];
        end
        [W1,Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(n1, n2, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability);
        
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

        % True if the Jaccard index passes a certain level, e.g. 0.8
        J = length(intersect(Clist_GT(1:cluster_size),Clist(1:CID)))/length(union(Clist_GT(1:cluster_size),Clist(1:CID)));
        if J >= thresh
            P_value = permutation_test_ttest(ctrl_mtx, case_mtx, Wp, CID, Clist, M_perm, threshold_GT, Method);
            power_list(m) = double(P_value<= FWER_threshold);
            %         disp(double(P_value<= FWER_threshold))
        else
            power_list(m) = 0;
        end
        send(dq, m);
    end
else
    % Serial execution with regular for loop
    for m = 1:M_rep
        thresh = 0.8;
        theta = Cohensd*sigma0;
        L = chol(Covariance);
        [W1,Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(n1, n2, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability);
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

        % True if the Jaccard index passes a certain level, e.g. 0.8
        J = length(intersect(Clist_GT(1:cluster_size),Clist(1:CID)))/length(union(Clist_GT(1:cluster_size),Clist(1:CID)));
        if J >= thresh
            P_value = permutation_test_ttest(ctrl_mtx, case_mtx, Wp, CID, Clist, M_perm,threshold_GT, Method);
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