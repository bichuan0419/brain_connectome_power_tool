function P_value = permutation_test_reg(Wm, X, Wp, CID, Clist, M, func,threshold_GT)


% edge_density = @(C) sum(sum(C))/(size(C,1)*(size(C,1)-1));
% node_density = @(C) sum(sum(C))/size(C,1);

T_orig =  custom_statistic(Wp, Clist, CID, threshold_GT);


T_vec = zeros(M,1);
for m=1:M
    
    X_perm = X(randperm(length(X)));
    [~, ptemp] = corr(X_perm, Wm);
    
    Wp_temp=squareform(-log10(ptemp));
    [~, Clist_temp, CID_temp] = greedy(Wp_temp, func,threshold_GT);
    [~,T_vec(m)] = custom_statistic(Wp_temp, Clist_temp, CID_temp, threshold_GT);
    
end

P_value = sum((T_vec - T_orig) >0 )/M;

end