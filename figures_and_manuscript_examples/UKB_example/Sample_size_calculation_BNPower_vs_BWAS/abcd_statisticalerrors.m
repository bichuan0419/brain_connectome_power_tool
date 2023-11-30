function [type1,type2,typem,types,types_pvals] = abcd_statisticalerrors(Corrs,Pvals,fscorrs,fspvals,iter)
% Generate Fig 2 of MP 
% Calculates error rates for a given brain-phenotype correlation
% type1 = 0;
% type2 = 0;
% iter;
%% Type 1 error rate % Should in theory be whatever your pvalue is 

thr = [.0000001];
type1 = zeros(size(Corrs,3),length(thr));
for p = 1:length(thr)
    EffectsMats = fspvals;
    pthr = thr(p); % Bonferroni correction 
    EffectsMats(EffectsMats>=pthr) = 1;
    EffectsMats(EffectsMats<pthr) = 0;
    disp(num2str(p))
    Type1error = nan(size(Corrs,1),16,length(thr));
    for r = 1:size(Corrs,1) % Loop thru edges
        if EffectsMats(r,1) == 1 % if not significant 
            for bin = 1:size(Corrs,3)
                thesep = Pvals(r,:,bin);
                Type1error(r,bin,p) = 100*(length(find(thesep<thr(p) ))/iter);
            end
        end
    end
    thisthr = squeeze(Type1error(:,:,p));
    thisthr = thisthr(:, ~all(isnan(thisthr)));
    thisthr(isnan(thisthr)) = 0;
    type1(:,p) = squeeze(mean(thisthr,1))';
end


%% Type 2 error rate 
 

Type2error = nan(size(Corrs,1),size(Corrs,3),length(thr));
type2 = zeros(size(Corrs,3),length(thr));
for p = 1:length(thr)
    EffectsMats = fspvals;
    pthr = thr(p);
    EffectsMats(EffectsMats>=pthr) = 1;
    EffectsMats(EffectsMats<pthr) = 0;
    for r = 1:size(Corrs,1) % Loop thru edges
        if EffectsMats(r,1) == 0 % if significant 
            for bin = 1:size(Corrs,3)
                thesep = Pvals(r,:,bin);
                %thesec = Corrs(r,:,bin);
                %thesec(thesec>0) = 1;
                %thesec(thesec<0) = -1; %& thesec == directionidx(r,:) 
                Type2error(r,bin,p) = 100*(length(find(thesep>thr(p) ))/iter);
            end
        end
    end
    thisthr = squeeze(Type2error(:,:,p));
    type2(:,p) = squeeze(nanmean(thisthr,1))';
end

%% Type M - inflated correlation 
pvalthr = [.0000001];

% Index sig corrs, for each iteration find % deviating by >|0.10|
AllInflatedCorrs = cell(length(pvalthr),1);
for p = 1:length(pvalthr)
    thispthr = pvalthr(p);
    EffectsMats = fspvals;
    EffectsMats(EffectsMats>=thispthr) = 1;
    EffectsMats(EffectsMats<thispthr) = 0;

    idx = find(EffectsMats==0);

    % determine rate of inflated correlation for each significant correlation 
    % run of iteration of what we consider "inflated" (20%:20:200)
    thr = [1.2:.1:3];
    InflatedCorrs = zeros(length(idx),size(Corrs,3),length(thr));
    for t = 1:length(thr)
        ip = thr(t); % inflation threshold (e.g., 2x)
        for i = 1:length(idx)
           thiscorr = fscorrs(idx(i));
           allstudycorrs = squeeze(Corrs(idx(i),:,:));
           allstudypvals = squeeze(Pvals(idx(i),:,:));
           for bin = 1:size(allstudycorrs,2)
               thisbinc = allstudycorrs(:,bin);
               thisbinp = allstudypvals(:,bin);
               % zero out significant but wrong direction (type 1 error)
               if thiscorr > 0 
                   thisbinc(thisbinc<0) = 0; 
                   thisbinc(thisbinc>0 & thisbinp > thispthr) = 0;
               else
                   thisbinc(thisbinc>0) = 0; 
                   thisbinc(thisbinc<0 & thisbinp > thispthr) = 0;
               end
               if thiscorr > 0  % true positive corelation 
                   InflatedCorrs(i,bin,t) = 100*(length(find(thisbinc > thiscorr * ip)) / length(find(thisbinc~=0)));
               else % true negative correlation 
                   InflatedCorrs(i,bin,t) = 100*(length(find(thisbinc < thiscorr * ip)) / length(find(thisbinc~=0)));
               end
           end
        end
    end
    AllInflatedCorrs{p} = InflatedCorrs;
end
typem = AllInflatedCorrs;


%% Type S error - Sign Flip  

% all correlations 
directionidx = zeros(length(fscorrs),1);
directionidx(fscorrs>0) = 1;
directionidx(fscorrs<0) = -1;

Signflip = zeros(length(directionidx),size(Corrs,3));
for t = 1:size(Signflip,1)
   allstudycorrs = squeeze(Corrs(t,:,:));
   for bin = 1:size(allstudycorrs,2)
       thisbinc = allstudycorrs(:,bin);
       if directionidx(t,1) == 1 % positive correlation 
           Signflip(t,bin) = 100*(length(find(thisbinc < 0 )) / length(thisbinc));
       else % negative correlation 
           Signflip(t,bin) = 100*(length(find(thisbinc > 0 )) / length(thisbinc));
       end
   end
end

% sort correlations for plotting 
[~,sortidx] = sort(abs(fscorrs),'descend');
types = Signflip(sortidx,:);

% Loop thru p thresholds
types_pvals = cell(length(pvalthr),1);
for p = 1:length(pvalthr)
    idx = find(fspvals<pvalthr(p));
    newfs = fscorrs(idx);
    directionidx = zeros(length(newfs),1);
    directionidx(newfs>0) = 1;
    directionidx(newfs<0) = -1;

    Signflip = zeros(length(directionidx),size(Corrs,3));
    for t = 1:size(Signflip,1)
       allstudycorrs = squeeze(Corrs(idx(t),:,:));
       allstudypvals = squeeze(Pvals(idx(t),:,:));
       for bin = 1:size(allstudycorrs,2)
           thisbinc = allstudycorrs(:,bin);
           thisbinp = allstudypvals(:,bin);
           if directionidx(t,1) == 1 % positive correlation 
               Signflip(t,bin) = 100*(length(find(thisbinc < 0 & thisbinp < pvalthr(p))) / length(thisbinc));
           else % negative correlation 
               Signflip(t,bin) = 100*(length(find(thisbinc > 0 & thisbinp < pvalthr(p))) / length(thisbinc));
           end
       end
    end

    % sort correlations for plotting 
    [~,sortidx] = sort(abs(newfs),'descend');
    types_pvals{p} = Signflip(sortidx,:);
end





end

