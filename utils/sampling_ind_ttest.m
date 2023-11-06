function [W1,Wp,Clist_GT, case_mtx, ctrl_mtx, threshold_GT] = sampling_ind_ttest(ctrl_num, case_num, N, sigma0, theta, cluster_size, rho_in, rho_out, L, Reliability)
            % a simple data generator
            % Chuan

            % output: W1 ground truth
            %         Wp permuted data

            % Load parameters
            Z = pick_idx(cluster_size, N);

            edge_num = N*(N-1)/2;
            all_idx = 1:edge_num;
            case_idx = Z;
            % non-covariate-related edges
            ctrl_idx = all_idx;
            ctrl_idx(case_idx) = [];
            % add non central dist. to the case
            % calculate the number of FP edges
            d_out = floor(rho_out * length(ctrl_idx));
            % number of FN edges
            d_in = floor((1-rho_in) * length(case_idx));
            % obtain the indices
            d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
            d_out_orig = ctrl_idx(d_out_idx);
            ctrl_idx(d_out_idx) = [];
            d_in_idx = randsample(1:length(case_idx),d_in,false);
            d_in_orig = case_idx(d_in_idx);
            case_idx(d_in_idx) = [];
            
            Reliability_vec = sqrt(squareform(Reliability));
            theta_vec = theta.* Reliability_vec;
            case_mtx = zeros(case_num, edge_num);
            case_mtx(:,case_idx) = randn(case_num, length(case_idx))*sigma0 + theta_vec(case_idx);
            case_mtx(:,d_out_orig) = randn(case_num, d_out)*sigma0 + theta_vec(d_out_orig);
            case_mtx(:,ctrl_idx) = randn(case_num, length(ctrl_idx))*sigma0;
            case_mtx(:,d_in_orig) = randn(case_num, d_in)*sigma0;
            ctrl_mtx = randn(ctrl_num, edge_num)*sigma0;
            
            case_mtx = case_mtx*L;
            ctrl_mtx = ctrl_mtx*L;
            
            pm = randperm(N);
            [~,p_vec]=ttest2(case_mtx,ctrl_mtx);
            W1=squareform(-log10(p_vec));
            Wp = W1(pm,pm);
            % return ground truth Clist
            [~,Clist_GT] = sort(pm);
            
            %% find a good cut
            threshold_ratio = 1:20;
            threshold_vec = prctile(p_vec,threshold_ratio);
            f1score_vec = zeros(length(threshold_vec),1);
            for i = 1:length(threshold_vec)
                target = ones(1,edge_num);
                target(ctrl_idx)=0;
                output = p_vec <= threshold_vec(i);
                cm = confusionmat(double(target),double(output));
                f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
            end

            [~,thresh_idx] = max(f1score_vec);
            threshold_GT=threshold_vec(thresh_idx);

        end