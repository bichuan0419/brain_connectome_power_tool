function [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, DF, cluster_size, rho_in, rho_out, L, Reliability)
% SAMPLING_IND_REG Simulates data for independent regression analysis.
%
% This function generates a simulated dataset for regression analysis,
% where the network structure, covariate-related effects, and non-covariate-related
% effects are taken into consideration. It randomly selects covariate-related edges
% within a network and assigns different distributions to covariate-related and
% non-covariate-related edges based on specified parameters.
%
% Inputs:
%   subject_num: Integer, the number of subjects in the study.
%   N: Integer, the size of the entire network (number of nodes).
%   sigma0: Double, standard deviation of the non-covariate-related edges.
%   f2: Double, Cohen's f2 effect size between the independent variable and
%       covariate-related edges.
%   DF: Double, Degrees of freedom.
%   cluster_size: Integer, size of the cluster within the network where the
%                 covariate-related effects are concentrated.
%   rho_in: Double, proportion of significant edges within the covariate-related
%           subnetworks.
%   rho_out: Double, proportion of significant edges within the non-covariate-related
%            edges.
%   L: Matrix, Cholesky decomposition of the desired correlation structure for
%      simulated data.
%   Reliability: Matrix, reliability coefficients for each edge.
%
% Outputs:
%   W1: Matrix, -log10 transformed p-value adjacency matrix for the original data.
%   Wp: Matrix, permuted -log10 p-value adjacency matrix.
%   Clist_GT: Vector, list of nodes representing the ground truth clustering.
%   Wm: Matrix, raw data matrix for the non-covariate-related edges.
%   X: Vector, sampled independent variable values for subjects.
%   threshold_GT: Double, determined threshold for the ground truth based on F1 score.
%
% Author:
%   Chuan - Original implementation for data sampling in network analysis.
%
% See also:
%   CORR, RANDN, RANDPERM, SQUAREFORM


% Initialize variables and constants
% Pick indices for covariate-related clusters within the network
Z = pick_idx(cluster_size, N);
% Calculate the total number of possible edges in an undirected graph without self-loops
edge_num = N*(N-1)/2;
% Create an array of all possible edge indices
all_idx = 1:edge_num;
% Assign covariate-related edges based on picked indices
case_idx = Z;
% Initialize control indices for non-covariate-related edges
ctrl_idx = all_idx;
% Exclude covariate-related edges from control indices
ctrl_idx(case_idx) = [];

% Add non-central distribution to the covariate-related cases
% Determine the number of false positive edges based on the specified proportion
d_out = floor(rho_out * length(ctrl_idx));
% Determine the number of false negative edges based on the specified proportion
d_in = floor((1-rho_in) * length(case_idx));
% Obtain random indices for false positives within control edges
d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
% Store the original control indices corresponding to false positives
d_out_orig = ctrl_idx(d_out_idx);
% Remove the selected false positive indices from control list
ctrl_idx(d_out_idx) = [];
% Obtain random indices for false negatives within case edges
d_in_idx = randsample(1:length(case_idx),d_in,false);
% Store the original case indices corresponding to false negatives
d_in_orig = case_idx(d_in_idx);
% Remove the selected false negative indices from case list
case_idx(d_in_idx) = [];
% Combine case indices with false positive indices to get polluted case indices
case_idx_polluted = [case_idx d_out_orig];
% Combine control indices with false negative indices to get polluted control indices
ctrl_idx_polluted = [ctrl_idx d_in_orig];

% Sampling distributions for covariate-related and non-covariate-related edges
% Generate independent variable X from a normal distribution with specified standard deviation
X = randn(subject_num, 1) * sigma0;
% Generate matrix for non-covariate-related edges with standard deviation sigma0
Wm = randn(subject_num, edge_num) * sigma0;
% Calculate reliability vector for the edges
Reliability_vec = sqrt(squareform(Reliability));
% Calculate f2 effect size vector adjusted by reliability for each edge
f2_vec = f2 .* Reliability_vec;
% Calculate the corresponding R^2 values for the f2 effect sizes
R2_vec = f2_vec ./ (1 + f2_vec);
% Adjust covariate-related edges based on R^2 values and independent variable X
% simulating the effect size with added noise
Wm(:, case_idx_polluted) = repmat(sqrt(R2_vec(case_idx_polluted)), size(Wm, 1), 1) .* X ...
    + repmat(sqrt(1 - R2_vec(case_idx_polluted)), size(Wm, 1), 1) .* randn(subject_num, length(case_idx_polluted)) * sigma0;
% Apply the Cholesky decomposition L to introduce the desired correlation structure to Wm
if ~isempty(L)
    Wm = Wm * L;
end

% Obtain p-values for the correlation between X and each edge
[corr_vec, ~] = corr(X, Wm);
% Calculate the t-values using the degrees of freedom
t_vec = corr_vec ./ sqrt((1 - corr_vec.^2) / DF);
% Calculate the p-values for a two-tailed test
p_vec = 2 * (1 - tcdf(abs(t_vec), DF));
% Find the smallest non-zero element
min_nonzero_p = min(p_vec(p_vec > 0));
% Replace zeros with the smallest non-zero value
p_vec(p_vec == 0) = min_nonzero_p;
% Randomly permute the nodes to simulate a null distribution
pm = randperm(N);
% Create a -log10 transformed p-value adjacency matrix from the correlation p-values
W1 = squareform(-log10(p_vec));
% Apply the permutation to the adjacency matrix to get the permuted matrix
Wp = W1(pm, pm);
% Return the ground truth clustering list, representing the true node order
[~, Clist_GT] = sort(pm);

% Find a good threshold for statistical significance
threshold_ratio = logspace(-1,1,50); % Define a set of percentile thresholds
% Determine the threshold values corresponding to the percentiles in p-value distribution
threshold_vec = prctile(p_vec, threshold_ratio);
% Initialize a vector to store F1 scores for each threshold
f1score_vec = zeros(length(threshold_vec), 1);
% Loop over each threshold to calculate the F1 score
for i = 1:length(threshold_vec)
    % Define the ground truth for edges: 1 for covariate-related, 0 for non-covariate-related
    target = ones(1, edge_num);
    target(ctrl_idx) = 0;
    % Define the predicted edges based on the current threshold
    output = p_vec <= threshold_vec(i);
    % Calculate the confusion matrix for the current threshold
    cm = confusionmat(double(target), double(output));
    % Calculate the F1 score from the confusion matrix
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

[~,thresh_idx] = max(f1score_vec);
threshold_GT=threshold_vec(thresh_idx);

end