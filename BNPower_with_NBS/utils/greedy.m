function [W_greedy, Clist, CID] = greedy(Wp, func, threshold_GT)
% a greedy version of the greedy peeling algorithm
% the algorithm finds a cluster at each time until the average value in the
% remainder adjacency matrix is below the mean of the input cluster
% Chuan

Wp(Wp<-log10(threshold_GT)) = 0;
N = size(Wp,1);
Recording_Matrix = zeros(N,2);
Recording_Matrix(1,:) = [0  func(Wp)];
Recording_Clist = 1:N;
Wp_temp = Wp;
iter = 1;
while iter < N
    idxlist_temp = 1:length(Wp_temp);
    % work on the temp matrix
    [~,idx_min_temp] = min(sum(Wp_temp));
    % find corresponding index for idx_min_temp in Recording_Clist
    idxlist_temp(idx_min_temp) = [];
    
    score_temp = func(Wp_temp(idxlist_temp,idxlist_temp));
    % record index and score
    Recording_Matrix(iter,:) = [Recording_Clist(idx_min_temp),score_temp];
    Recording_Clist(idx_min_temp) = [];
    % update matrix
    Wp_temp = Wp_temp(idxlist_temp,idxlist_temp);
    iter = iter+1;
end
Recording_Matrix(N,1) = Recording_Clist;
[~,max_idx] = max(Recording_Matrix(:,2));
if max_idx == 1
    Clist = Recording_Matrix(end:-1:1,1);
    W_greedy = Wp(Clist,Clist);
    CID = N;
else
    removing_node = Recording_Matrix(max_idx:-1:1,1);
    Node_Seq = Recording_Matrix(end:-1:max_idx+1,1);
    Clist = [Node_Seq;removing_node];
    W_greedy = Wp([Node_Seq;removing_node],[Node_Seq;removing_node]);
    CID = length(Node_Seq);
end

end