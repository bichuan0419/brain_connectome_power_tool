function P_value = permutation_test_ttest_NBS(ctrl_mtx, case_mtx, CID, M, threshold_GT)
    % Conducts a permutation test using the t-test to assess the statistical significance
    % of a test statistic derived from control and case matrices compared to a permuted dataset.
    
    % Inputs:
    % ctrl_mtx    - Matrix containing control group data.
    % case_mtx    - Matrix containing case group data.
    % Wp          - A matrix similar to the input matrices but with permuted edges.
    % CID         - Cluster IDs for calculating the custom statistic.
    % Clist       - A list of clusters corresponding to the cluster IDs.
    % M           - The number of permutations to be carried out.
    % func        - A function handle for the clustering algorithm used in the permutation test.
    % threshold_GT- The threshold used for defining significant edges in the test statistic.

    % Outputs:
    % P_value     - The p-value representing the significance level of the observed test statistic.
    
    %% FOR NBS, the test statistics is changed to the size of the largest subnetwork
    % Calculate the original test statistic using the weighted permuted matrix (Wp), 
    % the list of clusters (Clist), cluster IDs (CID), and the ground truth threshold.
    T_orig = max(CID);

    % Initialize an array to hold the test statistics from each permutation.
    T_vec = zeros(M, 1);
    
    % Concatenate control and case matrices for permutation.
    Wt_orig = [ctrl_mtx; case_mtx];

    % Start permutation loop.
    for m = 1:M
        % Permute the rows of the concatenated matrix.
        Wt_c = Wt_orig(randperm(size(Wt_orig, 1)), :);
        % Separate the permuted matrix into new control and case matrices.
        Wt1_c = Wt_c(1:size(ctrl_mtx, 1), :);
        Wt2_c = Wt_c((size(ctrl_mtx, 1) + 1):end, :);
        % Perform the t-test on the permuted control and case matrices.
        [~, ptemp] = ttest2(Wt1_c, Wt2_c);
        
        % Convert the p-values to a matrix format and apply negative log transformation.
        Wp_temp = squareform(-log10(ptemp));
        % Use the NBS algorithm to find clusters in the permuted matrix.
        [~, ~, CID_temp] = NBS(Wp_temp, -log10(threshold_GT));
        % Calculate the test statistic for this permutation.
        T_vec(m) = max(CID_temp);
    end

    % Determine the p-value as the proportion of permuted statistics that exceed
    % the original test statistic, indicating significance.
    P_value = sum((T_vec - T_orig) > 0) / M;

end
