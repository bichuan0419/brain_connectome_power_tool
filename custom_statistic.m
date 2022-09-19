function [T0_vec, score_max] = custom_statistic(Wp, Clist, CID, threshold_GT)
            % calculate test statistics from permutation test
            
            all_pvec  = squareform(Wp);
            % get number of clusters
            cluster_len = length(CID);
            % define test statistic for each cluster
            T0_vec = zeros(cluster_len,1);
            % loop each cluster and calculate the statistic for each cluster
            for i = 1:cluster_len
                % get the current cluster of interest
                if i == 1
                    W_cluster = Wp(Clist(1:CID(1)),Clist(1:CID(1)));
                else
                    W_cluster = Wp(Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)), ...
                        Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)));
                end
                % neglect small clusters
                if size(W_cluster,1) <= 5
                    T0_vec(i) = 0;

                else
                    % initialize score for the current cluster
                    % calculate the statistic
                    if W_cluster == 0
                        score = 0;
                    else
                        cur_pvec = squareform(W_cluster);
                        cur_pvec(cur_pvec < -log(threshold_GT)) = 0;
                        cur_pvec(cur_pvec >= -log(threshold_GT)) = 1;
%                         score = sum(cur_pvec)/(length(W_cluster) * (length(W_cluster)-1));
                    score = (sum(cur_pvec)/(length(W_cluster) * (length(W_cluster)-1))) * (sum(cur_pvec)/length(all_pvec(all_pvec > threshold_GT)));
                    end
                    T0_vec(i) = score;

                end

            end
            % get the maximum of the scores across all clusters
            [score_max,id_max] = max(T0_vec);
        end
