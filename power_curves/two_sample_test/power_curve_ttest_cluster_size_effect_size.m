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
N = 90;
%% cluster_size_list
% proportion of significant edges of the covariate-related subnetwork
cluster_size_list = [5, 10, 15, 20, 25, 30, 40, 50];
% difference between the mean values mu1 - mu0  of the two distributions
theta_list = 0.2:0.1:1;
power_ttest_cluster_size = zeros(length(theta_list), length(cluster_size_list));
rho_out = 0.02;
rho_in = 0.5;
for i = 1:length(theta_list)
    for j = 1:length(cluster_size_list)
        power_ttest_cluster_size(i,j) = covariate_net_ttest(n, n, N, sigma0, theta_list(i),...
                    rho_in, rho_out, cluster_size_list(j), FWER_threshold, M_rep, M_perm);
    end
end

%% spline curve fitting
xq_cohensd = linspace(min(theta_list), max(theta_list), 100);
yq_cohensd = zeros(100, length(cluster_size_list));
for j = 1:length(cluster_size_list)
    pp_ttest = pchip(theta_list, power_ttest_cluster_size(:,j));
    yq_cohensd(:,j) = ppval(pp_ttest, xq_cohensd);
end
xq_used = xq_cohensd;
yq_used = yq_cohensd(:,2:end);
figure;
% plot(xq_cohensd(xq_cohensd<0.7), yq_cohensd(xq_cohensd<0.7,3:end), 'LineWidth', 2)
plot(xq_used, yq_used, 'LineWidth', 2)
xlabel('Cohen''s d', 'FontSize',24,'FontWeight','bold');
ylabel('Power', 'FontSize',24,'FontWeight','bold');
hold on
% Add the horizontal line
line([xq_used(1), xq_used(end)], [0.8 0.8], 'Color', 'red', 'LineWidth', 4 , 'LineStyle', '--', 'DisplayName', '');
hold off
legend({['|V_c| = ', num2str(cluster_size_list(2))], ...
    ['|V_c| = ', num2str(cluster_size_list(3))], ['|V_c| = ', num2str(cluster_size_list(4))], ...
    ['|V_c| = ', num2str(cluster_size_list(5))], ['|V_c| = ', num2str(cluster_size_list(6))], ...
    ['|V_c| = ', num2str(cluster_size_list(7))], ['|V_c| = ', num2str(cluster_size_list(8))]}, 'FontSize',12,'Location','southeast')
set(gcf, 'position', [1000         858         579         480]);


