function [W_DSD_greedy, Clist,Node_Seq,  removing_node] = greedy_peeling(Wp_DSD, lambda)
% a simple implementation of the greedy algorithm from the Denser than
% dense paper, using the generalized objective function (SICERS
% Biostatistics by Chuo Chen, Chuan Bi, et al)
% note that this only extract ONE dense subgraphs
% Chuan

% Wp_DSD(Wp_DSD<threshold_DSD) = 0;
N = size(Wp_DSD,1);
Recording_Matrix = zeros(N,2);
Recording_Clist = 1:N;
Wp_temp = Wp_DSD;
obj_func = @(C) (sum(squareform(C))/length(squareform(C))).^lambda.*(sum(squareform(C))).^(1-lambda);

for i = N:-1:1
    idxlist_temp = 1:length(Wp_temp);
    % work on the temp matrix
    [~,idx_min_temp] = min(sum(Wp_temp));
    % find corresponding index for idx_min_temp in Recording_Clist
    idxlist_temp(idx_min_temp) = [];
    
    score_temp = obj_func(Wp_temp(idxlist_temp,idxlist_temp));
    % record index and score
    try
        Recording_Matrix(N-i+1,:) = [Recording_Clist(idx_min_temp),score_temp];
        Recording_Clist(idx_min_temp) = [];
        % update matrix
        Wp_temp = Wp_temp(idxlist_temp,idxlist_temp);
    catch
        Recording_Matrix(N-i+1,:) = [Recording_Clist(idx_min_temp),0];
        Recording_Clist(idx_min_temp) = [];
        % update matrix
        Wp_temp = Wp_temp(idxlist_temp,idxlist_temp);
    end
end
[~,max_idx] = max(Recording_Matrix(:,2));
removing_node = Recording_Matrix(1:max_idx,1);
Node_Seq = Recording_Matrix(end:-1:max_idx+1,1);
Clist = [Node_Seq;removing_node];
W_DSD_greedy = Wp_DSD(Clist,Clist);

end