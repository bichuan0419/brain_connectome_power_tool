function [T0_vec, score_max] = custom_statistic(Wp, Clist, CID, r_vec, weights)
            % calculate test statistics from permutation test
            
            weights = reshape(weights,1, length(weights));
            % get number of clusters
            cluster_len = length(CID);
            % define test statistic for each cluster
            T0_vec = zeros(cluster_len,1);
            % loop each cluster and calculate the statistic for each cluster
            for i = 1:cluster_len

                if i == 1
                    W_cluster = Wp(Clist(1:CID(1)),Clist(1:CID(1)));
                else
                    W_cluster = Wp(Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)), ...
                        Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)));
                end

                if size(W_cluster,1) <= 5
                    T0_vec(i) = 0;

                else
                    % initialize score for the current cluster
                    score = zeros(length(r_vec),1);
                    % calculate the statistic
                    for j = 1:length(r_vec)
                        if W_cluster == 0
                            score(j) = 0;
                        else
                            cur_pvec = squareform(W_cluster);
                            cur_pvec(cur_pvec < -log(r_vec(j))) = 0;
                            cur_pvec(cur_pvec >= -log(r_vec(j))) = 1;
                            score(j) = sum(cur_pvec)/(length(W_cluster));
                        end

                    end
                    T0_vec(i) = weights*score;

                end

            end
            % get the maximum of the scores across all clusters
            [score_max,id_max] = max(T0_vec);
        end
