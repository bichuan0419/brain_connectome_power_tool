%% power curve
clc;clear;close all;
%%
% This code plots power curve
% T-test, vary sample sizes
% Number subjects
n = 50;
% SD for normal distribution that governs the correlations of
% non-covariate-related edges
sigma0 = 1;
FWER_threshold = 0.05;
% Number of repetitions
M_rep = 100;
% Number of permutation tests in one repetition
M_perm = 100;
% Size of the entire network
N = 100;
%% rho_in_list
% proportion of significant edges of the covariate-related subnetwork
rho_in_list = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
% difference between the mean values mu1 - mu0  of the two distributions
f2_list = 0.01:0.03:0.3;
power_reg_rho_in = zeros(length(f2_list), length(rho_in_list));
rho_out = 0.02;
cluster_size = 25;
for i = 1:length(f2_list)
    for j = 1:length(rho_in_list)
        power_reg_rho_in(i,j) = covariate_net_reg(n, N, sigma0, f2_list(i), rho_in_list(j), rho_out, cluster_size, FWER_threshold, M_rep, M_perm);
    end
end

%% spline curve fitting
xq_cohensd = linspace(min(f2_list), max(f2_list), 100);
yq_cohensd = zeros(100, length(rho_in_list));
for j = 1:length(rho_in_list)
    pp_ttest = pchip(f2_list, power_reg_rho_in(:,j));
    yq_cohensd(:,j) = ppval(pp_ttest, xq_cohensd);
end
xq_used = xq_cohensd;
yq_used = yq_cohensd(:,2:end);

figure;
plot(xq_used, yq_used, 'LineWidth', 2)
% plot(xq_cohensd, yq_cohensd, 'LineWidth', 2)
xlabel('Cohen''s f^2', 'FontSize',24,'FontWeight','bold');
ylabel('Power', 'FontSize',24,'FontWeight','bold');
hold on
% Add the horizontal line
line([xq_used(1), xq_used(end)], [0.8 0.8], 'Color', 'red', 'LineWidth', 4 , 'LineStyle', '--', 'DisplayName', '');
hold off
legend({['|\rho_1| = ', num2str(rho_in_list(2))], ['|\rho_1| = ', num2str(rho_in_list(3))], ['|\rho_1| = ', num2str(rho_in_list(4))], ...
    ['|\rho_1| = ', num2str(rho_in_list(5))], ['|\rho_1| = ', num2str(rho_in_list(6))], ...
    ['|\rho_1| = ', num2str(rho_in_list(7))], ['|\rho_1| = ', num2str(rho_in_list(8))]}, 'FontSize',16,'Location','southeast')
set(gcf, 'position', [1000         858         579         480]);

%%
% proportion of significant edges of the outside the covariate-related subnetwork
rho_in_list = [0.05 0.04 0.03 0.02 0.01];
% number of nodes in the covariate-related subnetworks
cluster_size_list = [10 20 30 40 50];
