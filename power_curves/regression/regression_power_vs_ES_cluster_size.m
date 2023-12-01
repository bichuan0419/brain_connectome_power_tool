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
% Degrees of freedom
DF = 98;
% Covariance matrix
Covariance = [];
% Reliability matrix
Reliability = ones(N);
Reliability = Reliability - diag(diag(Reliability));
% Use parallel
Useparallel = 0;
%% rho_in_list
% proportion of significant edges of the covariate-related subnetwork
rho_in = 0.5;
% difference between the mean values mu1 - mu0  of the two distributions
f2_list = 0.01:0.03:0.3;
cluster_size_list = [5, 10, 15, 20, 25, 30, 40, 50];
power_reg_cluster_size = zeros(length(f2_list), length(cluster_size_list));
rho_out = 0.02;
for i = 1:length(f2_list)
    parfor j = 1:length(cluster_size_list)
        power_reg_cluster_size(i,j) = covariate_net_reg(n, N, sigma0,...
    f2_list(i), DF, rho_in, rho_out, cluster_size_list(j), Covariance, ...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel);
    end
end

%% spline curve fitting
xq_cohensf = linspace(min(f2_list), max(f2_list), 100);
yq_cohensd = zeros(100, length(cluster_size_list));
for j = 1:length(cluster_size_list)
    pp_reg = pchip(f2_list, power_reg_cluster_size(:,j));
    yq_cohensd(:,j) = ppval(pp_reg, xq_cohensf);
end
xq_used = xq_cohensf;
yq_used = yq_cohensd(:,3:end);

figure;
plot(xq_used, yq_used, 'LineWidth', 2)
% plot(xq_cohensd, yq_cohensd, 'LineWidth', 2)
xlabel('Cohen''s f^2', 'FontSize',24,'FontWeight','bold');
ylabel('Power', 'FontSize',24,'FontWeight','bold');
hold on
% Add the horizontal line
line([xq_used(1), xq_used(end)], [0.8 0.8], 'Color', 'red', 'LineWidth', 4 , 'LineStyle', '--', 'DisplayName', '');
hold off
legend({['|V_c| = ', num2str(cluster_size_list(3))], ['|V_c| = ', num2str(cluster_size_list(4))], ['|V_c| = ', num2str(cluster_size_list(5))], ...
    ['|V_c| = ', num2str(cluster_size_list(6))], ['|V_c| = ', num2str(cluster_size_list(7))], ...
    ['|V_c| = ', num2str(cluster_size_list(8))]}, 'FontSize',16,'Location','southeast')
set(gcf, 'position', [1000         858         579         480]);
% save 
save('regression_power_vs_ES_cluster_size.mat')
% save figure
saveas(gcf,'regression_power_vs_ES_cluster_size.eps','epsc');