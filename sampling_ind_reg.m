function [W1,Wp,Clist_GT,Wm, X,  threshold_GT] = sampling_ind_reg(subject_num, N, sigma0, f2, cluster_size, rho_in, rho_out)
% a simple data generator
% Chuan

% input:
% subject_num: number of subjects
% N: size of the entire network (nodes)
% sigma0: SD of the non-covariate-related edges
% f2: Cohen's f between the independent variable and
% covariate-related edges
% cluster_size: size of the cluster
% rho_in: proportion of significant edges within the
% covariate-related subnetworks
% rho_out: proportion of significant edges within the
% non-covariate-related edges

% output:
% X: sampled independent variable
% Wp: permuted -logp adjacency matrix

% randomly pick covariate-related edges and
% non-covariate-related edges
Z = pick_idx(cluster_size, N);
edge_num = N*(N-1)/2;
all_idx = 1:edge_num;
case_idx = Z;
% non-covariate-related edges
ctrl_idx = all_idx;
ctrl_idx(case_idx) = [];
% add non central dist. to the case
% calculate the number of FP edges
d_out = floor(rho_out * length(ctrl_idx));
% number of FN edges
d_in = floor((1-rho_in) * length(case_idx));
% obtain the indices
d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
d_out_orig = ctrl_idx(d_out_idx);
ctrl_idx(d_out_idx) = [];
d_in_idx = randsample(1:length(case_idx),d_in,false);
d_in_orig = case_idx(d_in_idx);
case_idx(d_in_idx) = [];
case_idx_polluted = [case_idx d_out_orig];
ctrl_idx_polluted = [ctrl_idx d_in_orig];

% obtain distributions for covariate-related distributions
% first we get an X
X = randn(subject_num, 1)*sigma0;
Wm = randn(subject_num, edge_num)*sigma0;  % non-covariate-related edges
Wm(:, case_idx_polluted) = sqrt(f2/(1+f2))*X + sqrt(1/(1+f2))*Wm(:, case_idx_polluted);  % covariate-related edges, r= sqrt(f2/(1+f2))

%% get p values
[corr_vec, p_vec] = corr(X, Wm);
pm = randperm(N);
W1=squareform(-log(p_vec));
Wp = W1(pm,pm);
% return ground truth Clist
[~,Clist_GT] = sort(pm);

%% find a good cut
threshold_vec = logspace(-5,-1,100);
f1score_vec = zeros(length(threshold_vec),1);
for i = 1:length(threshold_vec)
    target = ones(1,edge_num);
    target(ctrl_idx)=0;
    output = p_vec <= threshold_vec(i);
    cm = confusionmat(double(target),double(output));
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

[~,thresh_idx] = max(f1score_vec);
threshold_GT=threshold_vec(thresh_idx);

end