function power = covariate_net_test(n1, n2, N, sigma0, sigma1, theta, rho_in, rho_out, cluster_size, FWER_threshold, M_rep, M_perm)


func = @(C) sum(sum(C))/length(C)^1;
power_list = zeros(M_rep,1);
parfor m = 1:M_rep
    
    [W1,Wp, Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind(n1, n2, N, sigma0, sigma1, theta, cluster_size, rho_in, rho_out);
    % greedy
    [W_greedy, Clist_greedy, CID_greedy] = greedy(Wp, func,threshold_GT);
    % filter by Jaccard index
    J = length(intersect(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)))/length(union(Clist_GT(1:cluster_size),Clist_greedy(1:CID_greedy)));
    if J >= 0.5
        P_value = permutation_test(ctrl_mtx, case_mtx, Wp, CID_greedy, Clist_greedy, M_perm, func,threshold_GT);
        power_list(m) = double(P_value<= FWER_threshold);
%         disp(double(P_value<= FWER_threshold))
    else
        power_list(m) = 0;
    end
    %        figure;
    %     subplot(1,2,1)
    %     imagesc(W1);colorbar;colormap(jet);
    %     subplot(1,2,2)
    %     imagesc(W_greedy);colorbar;colormap(jet);
    %
    
end
power = mean(power_list);
end