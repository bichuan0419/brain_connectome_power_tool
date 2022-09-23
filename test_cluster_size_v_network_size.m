clc;clear;close all;
% Number of HC
n1 = 100;
% Number of Patients
n2 = 100;
% SD for normal distribution that governs the correlations of
% non-covariate-related edges
sigma0 = 1;
% SD for normal distribution that governs the correlations of
% Covariate-related edges
sigma1 =1;
% difference between the mean values mu1 - mu0  of the two distributions
theta =0.25;
% proportion of significant edges of the covariate-related subnetwork
rho_in = 0.7;
% proportion of significant edges of the outside the covariate-related subnetwork
rho_out = 0.01;
% number of nodes in the covariate-related subnetworks
cluster_size = 30;
FWER_threshold = 0.05;
% Number of repetitions
M_rep = 200;
% Number of permutation tests in one repetition
M_perm = 100;
% Size of the entire network
N = 100;

cluster_N_network = [18 60; 24 80; 30 100; 36 120; 42 140; 48 160; 54 180; 60 200; 66 220; 72 240];
elapsed_list = zeros(size(cluster_N_network,1),1);
power_list = zeros(size(cluster_N_network,1),1);
for i = 1:length(power_list)
[power, elapsed] = covariate_net_test(n1, n2, cluster_N_network(i,2), sigma0, sigma1, theta, rho_in, rho_out, cluster_N_network(i,1), FWER_threshold, M_rep, M_perm);
power_list(i) = power;
elapsed_list(i) = elapsed;
end

%% run time |V|^2
pp = polyfit(cluster_N_network(:,2), elapsed_list,2);
y1 = polyval(pp,cluster_N_network(:,2));
figure
plot(cluster_N_network(:,2),elapsed_list,'o')
hold on
plot(cluster_N_network(:,2), y1 )
hold off