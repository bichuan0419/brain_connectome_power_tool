function [W1,Wp,Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(ctrl_num, case_num, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability)
% sampling_ind_ttest - Simulates independent t-tests on data with specified reliability and effects
%
% Syntax:  [W1, Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(ctrl_num, case_num, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability)
%
% Inputs:
%    ctrl_num - Number of control samples
%    case_num - Number of case samples
%    N - Number of nodes in the network
%    sigma0 - Standard deviation of the underlying normal distribution
%    theta - Mean shift in the case group, affected by reliability
%    cluster_size - Number of nodes in the covariate-related cluster
%    rho_in - Proportion of false negatives (FN) in the case group
%    rho_out - Proportion of false positives (FP) in the control group
%    L - Cholesky decomposition matrix for inducing correlation
%    Reliability - Reliability of measurements in the network
%
% Outputs:
%    W1 - Matrix of -log10 transformed p-values from the true data
%    Wp - Matrix of -log10 transformed p-values from the permuted data
%    Clist_GT - Ground truth clustering list
%    case_mtx - Matrix of simulated case data
%    ctrl_mtx - Matrix of simulated control data
%    threshold_GT - Threshold for significance determined by ground truth data

% Author:
%   Chuan - Original implementation for data sampling in network analysis.


% Load parameters
Z = pick_idx(cluster_size, N); % Picks indices for a cluster of given size within N items

% Calculate the total number of edges in a fully connected graph of N nodes
edge_num = N*(N-1)/2;
all_idx = 1:edge_num; % Vector of all possible edge indices
case_idx = Z; % Indices of edges that are related to the case condition

% Identify control indices by excluding case indices from all indices
ctrl_idx = all_idx;
ctrl_idx(case_idx) = []; % Remove case indices to leave only control indices

% Determine the number of false positive edges based on specified rate (rho_out)
d_out = floor(rho_out * length(ctrl_idx)); 

% Determine the number of false negative edges based on specified rate (rho_in)
d_in = floor((1-rho_in) * length(case_idx));

% Randomly select indices to simulate false positives
d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
d_out_orig = ctrl_idx(d_out_idx); % Store original indices of false positives
ctrl_idx(d_out_idx) = []; % Remove false positive indices from control indices

% Randomly select indices to simulate false negatives
d_in_idx = randsample(1:length(case_idx),d_in,false);
d_in_orig = case_idx(d_in_idx); % Store original indices of false negatives
case_idx(d_in_idx) = []; % Remove false negative indices from case indices

% Create reliability vector, which will modulate the effect size based on reliability
Reliability_vec = sqrt(squareform(Reliability)); % Convert reliability to a vector form
theta_vec = theta.* Reliability_vec; % Scale the mean shift by reliability

% Initialize matrices for case and control data
case_mtx = zeros(case_num, edge_num); % Matrix for case group data
ctrl_mtx = randn(ctrl_num, edge_num)*sigma0; % Matrix for control group data with normally distributed random values

% Add the effects to the case matrix based on the indices determined above
case_mtx(:,case_idx) = randn(case_num, length(case_idx))*sigma0 + theta_vec(case_idx); % Add true effect to case group
case_mtx(:,d_out_orig) = randn(case_num, d_out)*sigma0 + theta_vec(d_out_orig); % Add false positives to case group
case_mtx(:,ctrl_idx) = randn(case_num, length(ctrl_idx))*sigma0; % Fill remaining entries for control group
case_mtx(:,d_in_orig) = randn(case_num, d_in)*sigma0; % Add false negatives to case group

% Apply the correlation structure to both matrices
if ~isempty(L)
    case_mtx = case_mtx*L;
    ctrl_mtx = ctrl_mtx*L;
end

% Random permutation of N elements, used for data permutation
pm = randperm(N);

% Perform t-tests between case and control data and store the p-values
[~,p_vec] = ttest2(case_mtx,ctrl_mtx);

% Transform p-values to -log10 scale and organize them into a matrix
W1 = squareform(-log10(p_vec));

% Permute the matrix W1 according to the random permutation pm
Wp = W1(pm,pm);

% Sort the permutation to get the ground truth cluster list
[~,Clist_GT] = sort(pm);

%% Find an optimal threshold for distinguishing case vs. control edges
% Initialize vectors to store threshold values and F1 scores
threshold_ratio = logspace(-1,1,50); % Define a set of percentile thresholds
threshold_vec = prctile(p_vec,threshold_ratio); % Compute actual threshold values
f1score_vec = zeros(length(threshold_vec),1); % Vector to store F1 scores for each threshold

% Loop over each threshold value to calculate F1 score
for i = 1:length(threshold_vec)
    target = ones(1,edge_num); % Initialize target vector with ones for edges
    target(ctrl_idx) = 0; % Set control indices to zero in the target vector
    output = p_vec <= threshold_vec(i); % Determine which p-values fall below the current threshold
    cm = confusionmat(double(target),double(output)); % Compute confusion matrix for the current threshold
    % Calculate F1 score from the confusion matrix
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

% Find the index of the maximum F1 score and use it to set the optimal threshold
[~,thresh_idx] = max(f1score_vec);
threshold_GT = threshold_vec(thresh_idx); % The optimal threshold for distinguishing case vs. control edges


end