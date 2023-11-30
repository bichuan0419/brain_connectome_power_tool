function [W_NBS, Clist_NBS, CID] = NBS(Wp, threshhold_NBS)
% a simple implementation of the network-based statistics
% Chuan

Wp_NBS = Wp;
Wp_NBS(Wp_NBS<threshhold_NBS) = 0;
% Wp_NBS(Wp_NBS>=threshhold_NBS) = 1;
% find Connected graph components
[Cindx,comp_sizes] = conncomp(graph(Wp_NBS));

[CID, C_size] = sort(comp_sizes,'descend');
Clist_NBS = [];
for i = 1:length( C_size)
    ith_comp_idx_list = find(Cindx ==  C_size(i));
    Clist_NBS = [Clist_NBS ith_comp_idx_list];
end
W_NBS = Wp_NBS(Clist_NBS,Clist_NBS);

end