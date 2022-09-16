clc;clear;close all;
% Number of HC
n1 = 40;
% Number of Patients
n2 = 40;
% STD of HC normal distribution
sigma0 = 1;
% STD of Patients normal distribution
sigma1 =1;
% difference between the mean values mu1 - mu0
theta =0.5;
% proportion of significant edges of the covariate-related subnetwork
rho_in = 0.6;
% proportion of insignificant edges of the outside the covariate-related subnetwork
rho_out = 0.01;
% number of nodes in the covariate-related subnetworks
cluster_size = 30;
FWER_threshold = 0.05;
% Number of repetitions
M_rep = 100;
% Number of permutation tests in one repetition
M_perm = 100;
% Size of the entire network
N = 100;
% r_vec = [0.05 0.01 0.001 0.0005 0.0001];
func = @(C) sum(sum(C))/length(C)^1;
weights = [0 0 1 0 0];
weights = weights/sum(weights);
tic;
power_list = zeros(M_rep,1);
parfor m = 1:M_rep
    
    [W1,Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind(n1, n2, N, sigma0, sigma1, theta, cluster_size, rho_in, rho_out);
    % greedy
    [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func,threshold_GT);
    r_vec = [10 2 1 0.1 0.01] * threshold_GT;
    % filter by Jaccard index
    J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
    disp(J)
    if J >= 0.7
        P_value = permutation_test(ctrl_mtx, case_mtx, Wp, CID_greedy, Clist_greedy, M_perm, r_vec, weights, func,threshold_GT);
        power_list(m) = double(P_value<= FWER_threshold);
        disp(double(P_value<= FWER_threshold))
    else
        power_list(m) = 0;
    end
    %        figure;
    %     subplot(1,2,1)
    %     imagesc(W1);colorbar;colormap(jet);
    %     subplot(1,2,2)
    %     imagesc(W_greedy);colorbar;colormap(jet);
    %
    
end
elapsed = toc;
disp(['Power is: ', num2str(mean(power_list)), ', run time: ', num2str(elapsed)]);
