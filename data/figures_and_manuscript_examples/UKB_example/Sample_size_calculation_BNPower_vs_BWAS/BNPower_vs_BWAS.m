clc;clear;close all
%% BNPower
% Retrieve input values from app components

% Simulation Parameters
M_rep = 100;       % Number of repetitions for the simulation
M_perm = 100;      % Number of permutations for the permutation test
N = 246;           % Total number of nodes
NE = N*(N-1)/2;    % Number of edges (calculated from total nodes)
n1 = 80;           % Number of subjects in group A
n2 = 80;           % Number of subjects in group B
sigma0 = 0.2;      % Standard deviation (sigma)
Cohensd = 0.4;     % Effect size (Cohen's d)
rho_in = 0.77;     % Proportion of significant edges within the covariate-related subnetwork
rho_out = 0.03;    % Proportion of significant edges outside the covariate-related subnetwork
cluster_size = 21; % Cluster size
FWER_threshold = 0.05; % Familywise error rate threshold

% Initialize matrices
Covariance_matrix = []; % Set to empty to save memory
Reliability_matrix = ones(N); % Initialize reliability matrix with ones
Reliability_matrix = Reliability_matrix - diag(diag(Reliability_matrix)); % Remove diagonal elements (set to 0)
power_list = zeros(M_rep,1); % Initialize power list for storing results

% Main simulation loop
for m = 1:M_rep
    % Calculate theta based on Cohen's d and sigma
    theta = Cohensd*sigma0;
    % Cholesky decomposition of the covariance matrix
    L = chol(Covariance_matrix);
    % Sampling and independent t-test
    [W1,Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(n1, n2, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability_matrix);
    
    % Greedy algorithm for network analysis
    func = @(C) sum(sum(C))/size(C,1); % Function for the greedy algorithm
    [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func, threshold_GT);

    % Calculate Jaccard index for the greedy algorithm results
    J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
    % Check if Jaccard index is above the threshold (e.g., 0.8)
    if J >= 0.8
        % Perform permutation test
        P_value = permutation_test_ttest(ctrl_mtx, case_mtx, Wp, CID_greedy, Clist_greedy, M_perm, func,threshold_GT);
        % Update power list based on permutation test results
        power_list(m) = double(P_value <= FWER_threshold);
    else
        power_list(m) = 0;
    end
end

% Calculate and display the average power
fprintf('The power using BNPower is: %.2f, with sample sizes SA = %d, SB = %d.\n', [mean(power_list), n1, n2]);
%% BWAS
% Initialize parameters for the simulation
% Total number of nodes in the network
N = 246;
NE = N*(N-1)/2; % Number of edges in the network (calculated from total nodes)
n1 = 7000; % Number of subjects in group A
n2 = 7000; % Number of subjects in group B
sigma0 = 0.2; % Standard deviation (sigma)
Cohensd = 0.4; % Effect size (Cohen's d)
rho_in = 0.77; % Proportion of significant edges within the covariate-related subnetwork
rho_out = 0.03; % Proportion of significant edges outside the covariate-related subnetwork
cluster_size = 21; % Size of the cluster in the network
Covariance_matrix = []; % Initialize covariance matrix as empty to save memory
Reliability_matrix = ones(N); % Initialize reliability matrix with ones
Reliability_matrix = Reliability_matrix - diag(diag(Reliability_matrix)); % Remove diagonal elements (set to 0)
power_list_BWAS = zeros(M_rep,1); % Initialize power list for storing results

% Set up parameters for the statistical test
theta = Cohensd*sigma0;
L = chol(Covariance_matrix); % Cholesky decomposition of the covariance matrix

% Sampling and independent t-test from population
[~, ~, ~, case_mtx, ctrl_mtx, ~] = sampling_ind_ttest(n1, n2, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability_matrix);

% Combine case and control matrices and create a factor vector
rmats = [case_mtx; ctrl_mtx]; % Combined matrix of case and control data
factor = [ones(n1, 1); zeros(n2, 1)]; % Factor vector with 1 for case and 0 for control

% Calculate correlations for each edge
for x = 1:size(rmats,2)
    [FullSampleCorrs(x,1), FullSamplePvals(x,1)] = corr(rmats(:,x), factor);
end

% Set number of bins for the analysis
bins = [1000];

% Perform edge-wise correlation with iterative reliability
[Corrs, Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(rmats, factor, bins, 100);

% Calculate error rates (Type I, Type II, etc.)
[type1, type2, typem, types, types_pvals] = abcd_statisticalerrors(Corrs, Pvals, FullSampleCorrs, FullSamplePvals, 100);

% Calculate mean Type II error rate
mtype2 = type2;

% Display the power of the BWAS
fprintf('The power using BWAS is: %.2f, with sample size S = %d.\n', [mean(100 - mtype2), bins(1)]);