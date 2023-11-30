clc;clear;close all
%%
% Read the matrix from the CSV file; first column is t-statistic, second is p-value
p_t_result = readmatrix("t_p_old_vs_young.csv"); 
% Extract the third column from the matrix
p_vec = p_t_result(:,3);
%%
% Calculate the negative log10 of the p-values
nlogp_vec = -log10(p_vec);
Wp = squareform(nlogp_vec);

% Plot the matrix as an image with a colorbar and set the colormap to 'jet'
figure;imagesc(Wp);colorbar;colormap("jet")
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)
saveas(gcf,'Wp_orig.png','png');
%% 
% Perform parameter tuning on the matrix Wp using a greedy algorithm
[lambda_out, cut_out] = param_tuning(Wp, 'greedy');
% Apply the greedy algorithm to finalize the network with the tuned parameters
[W_greedy, Clist_greedy, CID_greedy] = greedy_final(Wp, cut_out, lambda_out);
%% 
% Plot the resulting greedy network as an image with a colorbar and 'jet' colormap
figure;imagesc(W_greedy);colorbar;colormap("jet")
hold on
rectangle('Position', [1, 1, CID_greedy(1), CID_greedy(1)], 'EdgeColor', 'r', 'LineWidth', 3);
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)
saveas(gcf,'Wp_dense.png','png');

%% get edges
cluster = Clist_greedy(1:CID_greedy(1));
A = zeros(246);
A(cluster, cluster) = 1;
A = A - diag(diag(A));
A_vec = squareform(A);
[~,idx_A] = find(A_vec);

% Calculate the proportion of edges within the cluster above the cut-off
rho_in = sum(nlogp_vec(idx_A) > cut_out)/length(idx_A);

idx_B = setdiff(1:30135, idx_A);
% Calculate the proportion of edges outside the cluster above the cut-off
rho_out = sum(nlogp_vec(idx_B) > cut_out)/length(idx_B);

writematrix(idx_A,'idx_A.csv')
writematrix(idx_B,'idx_B.csv')

