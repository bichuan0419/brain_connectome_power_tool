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
theta =0.3;
% proportion of significant edges of the covariate-related subnetwork
rho_in = 0.6;
% proportion of significant edges of the outside the covariate-related subnetwork
rho_out = 0.01;
% number of nodes in the covariate-related subnetworks
cluster_size = 30;
FWER_threshold = 0.05;
% Number of repetitions
M_rep = 100;
% Number of permutation tests in one repetition
M_perm = 100;
% Size of the entire network
N = 220;

power = covariate_net_test(n1, n2, N, sigma0, sigma1, theta, rho_in, rho_out, cluster_size, FWER_threshold, M_rep, M_perm);