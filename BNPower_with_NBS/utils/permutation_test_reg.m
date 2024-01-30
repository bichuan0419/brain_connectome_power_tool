function P_value = permutation_test_reg(Wm, X, Wp, CID, Clist, M, threshold_GT, Method)
    % This function performs a permutation test to evaluate the significance of
    % a statistic derived from regression between a matrix Wm and a vector X,
    % against permuted data. It uses a specified threshold to calculate the statistic.

    % Inputs:
    % Wm          - A weighted matrix from the original data.
    % X           - A vector of variables that could potentially be related to the edges in Wm.
    % Wp          - A matrix similar to Wm but with permuted edges.
    % CID         - Cluster IDs that are used in the custom statistic calculation.
    % Clist       - A list of clusters or ordering that corresponds to CID.
    % M           - The number of permutations to perform.
    % threshold_GT- The ground truth threshold for the test statistic calculation.
    % Method      - Method to be used for the data-driven network-level
    %               power analysis.

    % Outputs:
    % P_value     - The p-value indicating the significance of the observed statistic.

    % Calculate the original test statistic using the permuted matrix, the cluster list,
    % cluster IDs, and the ground truth threshold.
    T_orig = custom_statistic(Wp, Clist, CID, threshold_GT);

    % Initialize a vector to store the test statistics from each permutation.
    T_vec = zeros(M, 1);

    % Begin the permutation loop.
    for m = 1:M
        % Randomly permute the entries of vector X.
        X_perm = X(randperm(length(X)));
        % Compute the correlation between the permuted vector and the original weighted matrix.
        [~, ptemp] = corr(X_perm, Wm);
        
        % Transform the p-values from the correlation into a matrix form and take the negative log10.
        Wp_temp = squareform(-log10(ptemp));
        % Select method
        if strcmp(Method,'BNPower')
            % Use the greedy algorithm to find clusters in the permuted matrix.
            func = @(C) sum(sum(C))/size(C,1);
            [~, Clist_temp, CID_temp] = greedy(Wp_temp, func, threshold_GT);
        elseif strcmp(Method,'NBS')
            % NBS for network analysis
            [~, Clist_temp, CID_temp] = NBS(Wp_temp, -log10(threshold_GT));
        end
        % Calculate the test statistic for the permuted data.
        [~, T_vec(m)] = custom_statistic(Wp_temp, Clist_temp, CID_temp, threshold_GT);
    end

    % Calculate the p-value as the proportion of times the permuted test statistic
    % is greater than the original test statistic.
    P_value = sum((T_vec - T_orig) > 0) / M;

end
