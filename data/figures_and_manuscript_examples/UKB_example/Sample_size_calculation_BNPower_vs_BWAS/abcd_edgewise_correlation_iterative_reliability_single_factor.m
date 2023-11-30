function [Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(rmats,factor,binsize,iter)
%UNTITLED Summary of this function goes here
%   Input:  rmats = 2D matrix (subject x edge)
%           factor = behavioral item to correlate with RSFC
%                     Should be n x 1, where n = subject length and 
%                     f = behavioral variables
%           binsize = size of subject bins (intervals of x subjects)
%           iter = number of iterations at each bin size

%   Output: Corrs = subject bin x iteration. Each entry is an
%           RSFC/behavior correlation. 
%           Pvals = subject bin x iteration. Each entry is a pvalue
%           for the associated RSFC/behavior correlation in NetworkCorrs.


% Loop thru subjects in interval of bin size, iter times at each bin.
Corrs = zeros(size(rmats,2),iter,length(binsize));
Pvals = zeros(size(rmats,2),iter,length(binsize));


for i = 1:length(binsize)
    tempC = zeros(size(rmats,2),iter);
    tempP = zeros(size(rmats,2),iter);
    for j = 1:1:iter
%         if mod(j,10) == 0
%             disp(['Loop ' num2str(i) ' iteration ' num2str(j)])
%         end
        %generate a random vector of length subject
        idx = datasample(1:size(rmats,1),binsize(i));
        % Index out matrices and factor
        TheseMats = rmats(idx,:);
        TheseBehav = factor(idx,1);
        TheseMats(isnan(TheseBehav),:) = [];
        TheseBehav(isnan(TheseBehav)) = [];
        for edge = 1:size(TheseMats,2) % loop thru each edge
            [tempC(edge,j), tempP(edge,j)] = corr(TheseMats(:,edge),TheseBehav);
        end
    end 
    Corrs(:,:,i) = tempC;
    Pvals(:,:,i) = tempP;  
end


end

