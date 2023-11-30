function [W_DSD_greedy, Clist, Node_Seq, removing_node] = OQC_greedyA(Wp_DSD, func)
% OQC_GREEDYA implements a greedy algorithm to extract a dense subgraph.
% This function is based on the approach described in the "Denser than Dense" paper.
% It extracts one dense subgraph from the given weighted adjacency matrix.
%
% Inputs:
%   Wp_DSD - Weighted adjacency matrix of the graph.
%   func - Function handle used to evaluate the quality of a subgraph.
%
% Outputs:
%   W_DSD_greedy - Reordered weighted adjacency matrix based on the extracted subgraph.
%   Clist - List of nodes in the extracted subgraph and the remaining nodes.
%   Node_Seq - Sequence of nodes in the extracted dense subgraph.
%   removing_node - Nodes that were removed to form the dense subgraph.

% Initialize variables
N = size(Wp_DSD, 1); % Number of nodes in the graph
Recording_Matrix = []; % Matrix to record the score of each subgraph
Recording_Clist = 1:N; % List of all nodes
Wp_temp = Wp_DSD; % Temporary working copy of the adjacency matrix

% Iteratively remove nodes to find the densest subgraph
for i = N:-1:1
    idxlist_temp = 1:length(Wp_temp); % Temporary index list for the current subgraph

    % Find the node with the minimum sum of edge weights
    [~, idx_min_temp] = min(sum(abs(Wp_temp)));

    % Remove the selected node from the temporary index list
    idxlist_temp(idx_min_temp) = [];

    % Calculate the score of the current subgraph without the removed node
    score_temp = func(Wp_temp(idxlist_temp, idxlist_temp));

    % Record the removed node and the score of the remaining subgraph
    Recording_Matrix = [Recording_Matrix; [Recording_Clist(idx_min_temp), score_temp]];

    % Update the list of remaining nodes
    Recording_Clist(idx_min_temp) = [];

    % Update the temporary adjacency matrix
    Wp_temp = Wp_temp(idxlist_temp, idxlist_temp);
end

% Identify the densest subgraph
[~, max_idx] = max(Recording_Matrix(:, 2));
removing_node = Recording_Matrix(1:max_idx, 1);
Node_Seq = Recording_Matrix(end:-1:max_idx+1, 1);

% Construct the list of nodes in the dense subgraph and the remaining nodes
Clist = [Node_Seq; removing_node];

% Reorder the original weighted adjacency matrix based on the dense subgraph
W_DSD_greedy = Wp_DSD([Node_Seq; removing_node], [Node_Seq; removing_node]);
end
