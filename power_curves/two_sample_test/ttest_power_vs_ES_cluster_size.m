%% power curve
clc;clear;close all;
%%
% This code plots power curve
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
% Covariance matrix
Covariance = [];
% Reliability matrix
Reliability = ones(N);
Reliability = Reliability - diag(diag(Reliability));
% Use parallel
Useparallel = 0;
%% cluster_size_list
% proportion of significant edges of the covariate-related subnetwork
cluster_size_list = [5, 10, 15, 20, 25, 30, 40, 50];
% difference between the mean values mu1 - mu0  of the two distributions
theta_list = 0.2:0.1:0.8;
% Cohen's d
Cohensd_list = theta_list*sigma0;
power_ttest_cluster_size = zeros(length(theta_list), length(cluster_size_list));
rho_out = 0.02;
rho_in = 0.5;
for i = 1:length(Cohensd_list)
    parfor j = 1:length(cluster_size_list)
        power_ttest_cluster_size(i,j) = covariate_net_ttest(n, n, N, sigma0, ...
    Cohensd_list(i), rho_in, rho_out, cluster_size_list(j), Covariance,...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel);
    end
end

%% spline curve fitting
xq_cohensd = linspace(min(Cohensd_list), max(Cohensd_list), 100);
yq_cohensd = zeros(100, length(cluster_size_list));
for j = 1:length(cluster_size_list)
    pp_ttest = pchip(Cohensd_list, power_ttest_cluster_size(:,j));
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
% save 
save('ttest_power_vs_ES_cluster_size.mat')
% save figure
saveas(gcf,'ttest_power_vs_ES_cluster_size.eps','epsc');

