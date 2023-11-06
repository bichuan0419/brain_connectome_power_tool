function power = covariate_net_reg(subject_num, N, sigma0,...
    f2, rho_in, rho_out, cluster_size, Covariance, ...
    FWER_threshold, M_rep, M_perm, Reliability)
% Setup the Total Iterations and Waitbar:
wb = waitbar(0, 'Please wait...');

% Check for Parallel Computing Toolbox:
hasParallelToolbox = license('test', 'Parallel_Computing_Toolbox') && ~isempty(ver('parallel'));

if hasParallelToolbox
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
if hasParallelToolbox

    parfor m = 1:M_rep
        L = chol(Covariance);
        [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, cluster_size, rho_in, rho_out, L, Reliability);
        % greedy
        func = @(C) sum(sum(C))/length(C);
        % [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func);

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
        send(dq, m);


    end
else
    for m = 1:M_rep
        L = chol(Covariance);
        [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, cluster_size, rho_in, rho_out, L, Reliability);
        % greedy
        func = @(C) sum(sum(C))/length(C);
        % [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func);

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