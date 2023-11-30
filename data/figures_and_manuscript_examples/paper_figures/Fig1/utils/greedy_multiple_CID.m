function [W_DSD_greedy, Clist, CID] = greedy_multiple_CID(Wp, func, threshold_DSD)
% GREEDY_MULTIPLE_CID performs a greedy clustering algorithm on a weighted graph.
% This function applies a greedy peeling algorithm to identify clusters in the graph.
% The algorithm iteratively finds clusters until the average value in the
% remaining adjacency matrix is below a specified threshold.
%
% Inputs:
%   Wp - Weighted adjacency matrix of the graph.
%   func - Function handle used to evaluate the quality of a cluster.
%   threshold_DSD - Threshold for determining significant edges in the graph.
%
% Outputs:
%   W_DSD_greedy - Reordered weighted adjacency matrix based on the clustering.
%   Clist - List of nodes in each cluster.
%   CID - Number of nodes in each identified cluster.

% Initialize the weighted adjacency matrix, setting insignificant edges to zero
Wp_DSD = Wp;
Wp_DSD(abs(Wp_DSD) < -log(threshold_DSD)) = 0;

% Initialize variables for the clustering process
Node_list = {}; % Store all node lists from each cluster
Clist = []; % List of nodes in each cluster
CID = []; % Number of nodes in each cluster
nc = 1;
orig_Node = (1:length(Wp))'; % Original node indices

% Continue clustering until almost all nodes are clustered
while length(Clist) < length(Wp) - 1
    % Perform greedy peeling to obtain the next cluster
    [~, Clist_temp, Node_Seq, remaining_node] = OQC_greedyA(Wp_DSD, func);

    % Update the adjacency matrix to remove the nodes of the current cluster
    Wp_DSD = Wp_DSD(remaining_node, remaining_node);

    % Record the nodes in the current cluster and update the cluster lists
    Node_list{nc} = orig_Node(Clist_temp(1:length(Node_Seq)));
    Clist = [Clist; orig_Node(Clist_temp(1:length(Node_Seq)))];
    CID = [CID; length(Node_Seq)];

    % Update the original node indices and increment the cluster count
    orig_Node = orig_Node(remaining_node);
    nc = nc + 1;
end

% Include any remaining uninformative nodes in the last cluster
Clist = [Clist; orig_Node];
CID = [CID; length(remaining_node)];

% Reorder the original weighted adjacency matrix based on the clustering
W_DSD_greedy = Wp(Clist, Clist);
end
