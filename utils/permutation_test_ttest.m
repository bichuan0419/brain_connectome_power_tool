function P_value = permutation_test_ttest(ctrl_mtx, case_mtx, Wp, CID, Clist, M, func,threshold_GT)


            % edge_density = @(C) sum(sum(C))/(size(C,1)*(size(C,1)-1));
            % node_density = @(C) sum(sum(C))/size(C,1);

            % switch between functions
            casenum = size(case_mtx,1);
            ctrlnum = size(ctrl_mtx,1);
            % Wp=res_summary.Wp;
            Wt_orig = [ctrl_mtx;case_mtx];
            
            T_orig =  custom_statistic(Wp, Clist, CID, threshold_GT);


            T_vec = zeros(M,1);
            for m=1:M

                Wt_c=Wt_orig(randperm(casenum+ctrlnum),:);
                Wt1_c=Wt_c(1:ctrlnum,:);
                Wt2_c=Wt_c(ctrlnum+1:end,:);
                [~,ptemp]=ttest2(Wt1_c,Wt2_c);

                Wp_temp=squareform(-log10(ptemp));
                [~, Clist_temp, CID_temp] = greedy(Wp_temp, func,threshold_GT);
                [~,T_vec(m)] = custom_statistic(Wp_temp, Clist_temp, CID_temp, threshold_GT);

            end

            P_value = sum((T_vec - T_orig) >0 )/M;

        end