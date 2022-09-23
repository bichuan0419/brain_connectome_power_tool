clc;clear;close all
% Number of subjects
subject_num =220;
% SD for normal distribution that governs the correlations of
% non-covariate-related edges
sigma0 = 1;
% Cohen's f
f2 = 0.1;
% proportion of significant edges of the covariate-related subnetwork
rho_in = 0.45;
% proportion of significant edges of the outside the covariate-related subnetwork
rho_out = 0.02;
% number of nodes in the covariate-related subnetworks
cluster_size = 20;
FWER_threshold = 0.05;
% Number of repetitions
M_rep = 100;
% Number of permutation tests in one repetition
M_perm = 100;
% Size of the entire network
N = 100;
%% regression test
tic;
func = @(C) sum(sum(C))/length(C)^1;
power_list = zeros(M_rep,1);
parfor m = 1:M_rep
    
    [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, cluster_size, rho_in, rho_out);
    % greedy
    [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func,threshold_GT);
    %     imagesc(W1);colorbar;colormap(jet);
    % filter by Jaccard index
    J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
    if J >= 0.5
        P_value = permutation_test_reg(Wm, X, Wp, CID_greedy, Clist_greedy, M_rep, func,threshold_GT);
        power_list(m) = double(P_value<= FWER_threshold);
        %         disp(double(P_value<= FWER_threshold))
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
power = mean(power_list);
elapsed = toc;
