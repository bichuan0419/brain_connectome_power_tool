%% power curve
clc;clear;close all;
%%
% This code plots power curve
% T-test, vary sample sizes
% effect size
Cohensd = 0.5;
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
rho_in_list = 0.1:0.1:0.8;
% difference between the mean values mu1 - mu0  of the two distributions
sample_size_list = 10:10:100;
power_ttest_rho_in = zeros(length(sample_size_list), length(rho_in_list));
rho_out = 0.02;
cluster_size = 25;
for i = 1:length(sample_size_list)
    for j = 1:length(rho_in_list)
        power_ttest_rho_in(i,j) = covariate_net_ttest(sample_size_list(i), sample_size_list(i), N, sigma0, Cohensd,...
                    rho_in_list(j), rho_out, cluster_size, FWER_threshold, M_rep, M_perm);
    end
end

%% spline curve fitting
xq_cohensd = linspace(min(sample_size_list), max(sample_size_list), 100);
yq_cohensd = zeros(100, length(rho_in_list));
for j = 1:length(rho_in_list)
    pp_ttest = pchip(sample_size_list, power_ttest_rho_in(:,j));
    yq_cohensd(:,j) = ppval(pp_ttest, xq_cohensd);
end
xq_used = xq_cohensd;
yq_used = yq_cohensd(:,3:end);

figure;
plot(xq_used, yq_used, 'LineWidth', 2)
% plot(xq_cohensd, yq_cohensd, 'LineWidth', 2)
xlabel('Sample Size', 'FontSize',24,'FontWeight','bold');
ylabel('Power', 'FontSize',24,'FontWeight','bold');
hold on
% Add the horizontal line
line([xq_used(1), xq_used(end)], [0.8 0.8], 'Color', 'red', 'LineWidth', 4 , 'LineStyle', '--', 'DisplayName', '');
hold off
legend({['\rho_1 = ', num2str(rho_in_list(3))], ['\rho_1 = ', num2str(rho_in_list(4))], ...
    ['\rho_1 = ', num2str(rho_in_list(5))], ['\rho_1 = ', num2str(rho_in_list(6))], ...
    ['\rho_1 = ', num2str(rho_in_list(7))], ['\rho_1 = ', num2str(rho_in_list(8))]}, 'FontSize',16,'Location','southeast')
set(gcf, 'position', [1000         858         579         480]);
