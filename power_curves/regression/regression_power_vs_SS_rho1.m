%% power curve
clc;clear;close all;
%%
% effect size
Cohensf2 = 0.15;
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
rho_in_list = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
cluster_size = 25;
% difference between the mean values mu1 - mu0  of the two distributions
sample_size_list = 10:10:100;
power_reg_cluster_size = zeros(length(sample_size_list), length(rho_in_list));
rho_out = 0.02;

for i = 1:length(sample_size_list)
    parfor j = 1:length(rho_in_list)
        power_reg_cluster_size(i,j) = covariate_net_reg(sample_size_list(i), N, sigma0,...
    Cohensf2, DF, rho_in_list(j), rho_out, cluster_size, Covariance, ...
    FWER_threshold, M_rep, M_perm, Reliability, Useparallel);
    end
end



%% spline curve fitting
xq_cohensd = linspace(min(sample_size_list), max(sample_size_list), 100);
yq_cohensd = zeros(100, length(rho_in_list));
for j = 1:length(rho_in_list)
    pp_ttest = pchip(sample_size_list, power_reg_cluster_size(:,j));
    yq_cohensd(:,j) = ppval(pp_ttest, xq_cohensd);
end
xq_used = xq_cohensd;
yq_used = yq_cohensd(:,4:end);

figure;
plot(xq_used, yq_used, 'LineWidth', 2)
% plot(xq_cohensd, yq_cohensd, 'LineWidth', 2)
xlabel('Sample Size', 'FontSize',24,'FontWeight','bold');
ylabel('Power', 'FontSize',24,'FontWeight','bold');
hold on
% Add the horizontal line
line([xq_used(1), xq_used(end)], [0.8 0.8], 'Color', 'red', 'LineWidth', 4 , 'LineStyle', '--', 'DisplayName', '');
hold off
legend({['|\rho_1| = ', num2str(rho_in_list(4))], ['|\rho_1| = ', num2str(rho_in_list(5))], ...
    ['|\rho_1| = ', num2str(rho_in_list(6))], ['|\rho_1| = ', num2str(rho_in_list(7))], ...
    ['|\rho_1| = ', num2str(rho_in_list(8))]}, 'FontSize',16,'Location','southeast')
set(gcf, 'position', [1000         858         579         480]);
% save 
save('regression_power_vs_SS_rho1.mat')
% save figure
saveas(gcf,'regression_power_vs_SS_rho1.eps','epsc');