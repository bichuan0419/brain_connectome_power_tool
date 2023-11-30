function [W_DSD_greedy, Clist, CID] = greedy_final(Wp, threshold_DSD, lambda)
% a greedy version of the greedy peeling algorithm
% the algorithm finds a cluster at each time until the average value in the
% remainder adjacency matrix is below the mean of the input cluster
% Chuan

% initialize
Wp_DSD = Wp;
Wp_DSD(Wp_DSD<threshold_DSD) = 0;

% store all node lists from each cluster
Node_list = {}; nc = 1;
Clist = []; CID = [];
orig_Node = (1:length(Wp))';
% while loop 
while length(Clist) < length(Wp)-1
    % do greedy peeling to get the next cluster
    [~, Clist_temp,Node_Seq,  remaining_node] = greedy_peeling(Wp_DSD, lambda);
    % update adjacency matrix
    Wp_DSD = Wp_DSD(remaining_node, remaining_node);
    % update 
    p_rem = mean(mean(Wp_DSD));
    % record nodes
    Node_list{nc} = orig_Node(Clist_temp(1:length(Node_Seq)));
    Clist = [Clist; orig_Node(Clist_temp(1:length(Node_Seq)))];
    CID = [CID; length(Node_Seq)];
    % updates
    orig_Node = orig_Node(remaining_node);
    nc = nc + 1;
end
% uninformative nodes
Clist = [Clist; orig_Node];
CID = [CID; length(remaining_node)];
W_DSD_greedy = Wp(Clist, Clist);
end