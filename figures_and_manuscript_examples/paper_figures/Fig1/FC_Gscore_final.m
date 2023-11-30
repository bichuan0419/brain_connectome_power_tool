clc;clear;close all
%%
% Read data from a CSV file and remove the first row and column
df_reduced_UKB = readmatrix("data/beta_t_p.csv");
df_reduced_UKB = df_reduced_UKB(2:end, 2:end);
% df_reduced_UKB = df_reduced_UKB(:, 2:end);
% Compute the negative logarithm of the p-values (column 5 of the data)
nlogp_gscore_UKB = -log10(df_reduced_UKB(:,5));
% Identify indices with infinite values in nlogp_gscore_UKB
inf_idx_UKB = find(isinf(nlogp_gscore_UKB));
% Define the number of nodes and calculate the number of edges
nodeLen = 246;
edgeLen = nodeLen*(nodeLen-1)/2;
% Identify indices that are not infinite
noninf_idx = setdiff(1:edgeLen, inf_idx_UKB);
% Replace infinite values with the maximum of non-infinite values
nlogp_gscore_UKB(inf_idx_UKB) = max(nlogp_gscore_UKB(noninf_idx));
% Extract t-scores from the data
t_gscore_UKB = df_reduced_UKB(:,4);
% Convert vector to a square matrix form for visualization
Wp_gscore_UKB = squareform(nlogp_gscore_UKB);
Wt_gscore_UKB = squareform(t_gscore_UKB);


%% tstats
% Visualize the matrix of t-scores
figure;
imagesc(Wt_gscore_UKB);colorbar;colormap(jet);caxis([-15 15]);
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)

%% visualize both networks
% Visualize the matrix of negative log p-values
figure;
imagesc(Wp_gscore_UKB);colorbar;colormap(jet);caxis([0,40])
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)

%%
% Define parameters for network analysis
lambda_out_UKB = 1; cut_out_UKB = 1e-7;
% Define a function for network analysis
func_UKB = @(C) sum(sum(C))/size(C,1)^lambda_out_UKB;
% Perform greedy community detection
[W_greedy_UKB, Clist_greedy_UKB, CID_greedy_UKB] = greedy_multiple_CID(Wp_gscore_UKB, func_UKB, cut_out_UKB);

% Save the results to a MAT file
save('Fig1/results/UKB_dense_files.mat','Wp_gscore_UKB', 'Clist_greedy_UKB', 'CID_greedy_UKB');

%%
% Visualize the reordered matrix based on community detection
figure;
imagesc(Wp_gscore_UKB(Clist_greedy_UKB,Clist_greedy_UKB));colorbar;colormap(jet);caxis([0,40])
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)
% save to Fig1.eps
saveas(gcf,'Fig1b.eps','epsc');
greedy_edges = [squareform(W_greedy_UKB(1:CID_greedy_UKB(1),1:CID_greedy_UKB(1))),...
    squareform(W_greedy_UKB(CID_greedy_UKB(1)+1:sum(CID_greedy_UKB(1:2)),CID_greedy_UKB(1)+1:sum(CID_greedy_UKB(1:2))))] ;

edge_density_greedy = sum(greedy_edges >= -log(cut_out_UKB))/length(greedy_edges);

%%
% Calculate edge density and other metrics for different brain networks
% default mode network (DMN)
DMN = [3,5,6,11,13,14,23,33,34,35,41,42,43,44,51,52,79,80,81,83,84,87,88,95,121,122,141,144,153,154,175,176,179,181,187,188];
% Visual network (VN)
VN = [105,106,108,112,113,114,119,120,135,136,151,152,182,189,190,191,192,193,194,195,196,197,198,199,200,202,203,204,205,206,207,208,209,210];
% Sensorimotor (SN)
SN = [9,10,53,54,57,58,59,60,66,67,68,71,72,73,74,75,76,131,132,145,146,149,155,156,157,158,160,161,162,163,164,171,172];
SN = SN(randperm(length(SN)));
% Dorsal Attention
DAN = [7,8,25,26,30,55,56,63,64,85,86,91,92,97,98,107,125,126,127,128,129,130,133,134,139,140,143,150,159,201];
% Ventral Attention
VAN = [2,15,37,38,39,40,61,62,65,123,124,167,168,169,170,173,174,180,183,184,185,186];
% Limbic
LN = [27,45,47,48,49,50,69,70,77,78,89,90,93,94,96,101,102,103,104,109,110,111,115,116,117,118];
% FrontalParietal
FPN = [1,4,12,16,17,18,19,20,21,22,24,28,29,31,32,36,46,82,99,100,137,138,142,147,148,166];
p_DMN = squareform(Wp_gscore_UKB(DMN,DMN));
p_VN = squareform(Wp_gscore_UKB(VN, VN));
p_SN = squareform(Wp_gscore_UKB(SN, SN));
p_DAN = squareform(Wp_gscore_UKB(DAN, DAN));
p_VAN = squareform(Wp_gscore_UKB(VAN, VAN));
p_LN = squareform(Wp_gscore_UKB(LN, LN));
p_FPN = squareform(Wp_gscore_UKB(FPN, FPN));
%% visualize
% construct predefined global network
subnetworks_nodes = [DMN, DAN, FPN,   LN, SN,  VAN,VN];
rest_nodes = setdiff(1:246, subnetworks_nodes);
Clist_predefined = [subnetworks_nodes, rest_nodes];
figure;
imagesc(Wp_gscore_UKB(Clist_predefined, Clist_predefined)); colorbar; colormap(jet);caxis([0,40])
ax=gca;
ax.FontSize = 20;
xlabel('Brain Region','FontSize',24)
ylabel('Brain Region','FontSize',24)
% save to Fig1.eps
saveas(gcf,'Fig1a.eps','epsc');
%% calculate the ratio of significant edges within greedy subnetworks
total_sig_edges = sum(nlogp_gscore_UKB > -log(cut_out_UKB));
ratio_sig_edges_greedy = sum(greedy_edges >= -log(cut_out_UKB))/total_sig_edges;
%% do the same for other networks
% DMN
edge_density_DMN = sum(p_DMN >= -log(cut_out_UKB))/length(p_DMN);
ratio_sig_edges_DMN = sum(p_DMN >= -log(cut_out_UKB))/total_sig_edges;
% VN
edge_density_VN = sum(p_VN >= -log(cut_out_UKB))/length(p_VN);
ratio_sig_edges_VN = sum(p_VN >= -log(cut_out_UKB))/total_sig_edges;
% SN
edge_density_SN = sum(p_SN >= -log(cut_out_UKB))/length(p_SN);
ratio_sig_edges_SN = sum(p_SN >= -log(cut_out_UKB))/total_sig_edges;
% DAN
edge_density_DAN = sum(p_DAN >= -log(cut_out_UKB))/length(p_DAN);
ratio_sig_edges_DAN = sum(p_DAN >= -log(cut_out_UKB))/total_sig_edges;
% VAN
edge_density_VAN = sum(p_VAN >= -log(cut_out_UKB))/length(p_VAN);
ratio_sig_edges_VAN = sum(p_VAN >= -log(cut_out_UKB))/total_sig_edges;
% LN
edge_density_LN = sum(p_LN >= -log(cut_out_UKB))/length(p_LN);
ratio_sig_edges_LN = sum(p_LN >= -log(cut_out_UKB))/total_sig_edges;
% FPN
edge_density_FPN = sum(p_FPN >= -log(cut_out_UKB))/length(p_FPN);
ratio_sig_edges_FPN = sum(p_FPN >= -log(cut_out_UKB))/total_sig_edges;

%% plot comparison between the ratios of significant edges for each pre-defined networks
% first plot: edge density
X = categorical({'Visual','Sensorimotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default','Data-Driven'});
Y = [edge_density_VN, edge_density_SN, edge_density_DAN, edge_density_VAN, edge_density_LN, edge_density_FPN, edge_density_DMN, edge_density_greedy];

% Define the 8 colors
colors = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'};
rgb_colors = zeros(length(colors), 3);
for i = 1:length(colors)
    rgb_colors(i,:) = sscanf(colors{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
% Create a bar plot and assign colors to each bar
figure;
b = bar(X, Y,'FaceColor','flat');
for i = 1:8
    b.CData(i,:) = rgb_colors(i,:);
end
% Set fontsize for both x and y axis
set(gca,'FontSize',20);
ylabel('Sig. Edge Density','FontSize',20);
set(gcf, 'position',[546   415   476   599]);
% save to Fig1.eps
saveas(gcf,'Fig1c.eps','epsc');
%%
% second plot
X = categorical({'Visual','Sensorimotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default','Data-Driven'});
Y = [ratio_sig_edges_VN, ratio_sig_edges_SN, ratio_sig_edges_DAN, ratio_sig_edges_VAN, ratio_sig_edges_LN, ratio_sig_edges_FPN, ratio_sig_edges_DMN, ratio_sig_edges_greedy];

% Define the 8 colors
colors = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'};
rgb_colors = zeros(length(colors), 3);
for i = 1:length(colors)
    rgb_colors(i,:) = sscanf(colors{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
% Create a bar plot and assign colors to each bar
figure;
b = bar(X, Y,'FaceColor','flat');
for i = 1:8
    b.CData(i,:) = rgb_colors(i,:);
end
% Set fontsize for both x and y axis
set(gca,'FontSize',20);
ylabel('Identified Sig. Edges Ratio','FontSize',20);
set(gcf, 'position',[546   415   476   599]);
% save to Fig1.eps
saveas(gcf,'Fig1d.eps','epsc');