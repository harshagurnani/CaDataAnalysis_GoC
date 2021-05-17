function [pkcorr, zerocorr, corrs, lags, pkcorr2, zerocorr2, corrs2, pval2, event_indx,event_vals,event_times] = find_event_corrs( dFF, dFF_time, varargin )
    

DS_rate = [];
Perc_thresh=[];
Frac_thresh=[];
corr_bin = [];

if nargin > 2  

    charctr=0;
    numctr=0;
    start_read_var = 3;
    while start_read_var <= nargin
        vn = start_read_var-2;
       if ischar( varargin{vn} )
           charctr=1;
           if strcmp(varargin{vn}, 'ds_rate' )
               DS_rate = varargin{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargin{vn}, 'perc_thresh' )
               Perc_thresh = varargin{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargin{vn}, 'perc_thresh' )
                Frac_thresh = varargin{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargin{vn}, 'corr_bin' )
               corr_bin = varargin{vn+1};
               start_read_var = start_read_var + 2;
           else 
                error('Argument name %s is not an optional argument', varargin{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   DS_rate = varargin{vn};
                   start_read_var = start_read_var+1;
               case 2
                   Perc_thresh = varargin{vn};
                   start_read_var = start_read_var+1;
               case 3
                   Frac_thresh = varargin{vn};
                   start_read_var = start_read_var+1;
               case 4
                   corr_bin = varargin{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end

if isempty(DS_rate)
    DS_rate = 5;
end

if isempty(Perc_thresh)
    Perc_thresh = 98;
end


if isempty(Frac_thresh)
    Frac_thresh = 0.8;
end
if isempty(corr_bin)
    corr_bin = 150;
end


Tm = size(dFF,1);
nTrials = size(dFF,2);
nROI = size(dFF,3);

temp_dff = nan( ceil(Tm/DS_rate), nTrials, nROI);
temp_time = temp_dff;

parfor roi = 1:nROI
    for trial=1:nTrials
        temp_dff(:, trial, roi) = decimate( dFF(:,trial,roi), DS_rate);
        temp_time(:, trial, roi) = downsample(dFF_time(:,trial,roi), DS_rate);
    end
end

event_indx = cell(nROI, nTrials);
parfor roi = 1:nROI
   all_f = squeeze(temp_dff(:,:,roi));
   all_t = squeeze(temp_time(:,:,roi));
   pthresh = Frac_thresh * prctile( all_f(2:end)-all_f(1:end-1), Perc_thresh);

   for trial = 1:nTrials
        diffs = all_f(2:end,trial) - all_f(1:end-1,trial);
        event_indx{roi}{trial} = 1+find(diffs>=pthresh);
        event_times{roi}{trial} = all_t(event_indx{roi}{trial}, trial ); 
        event_vals{roi}{trial} = all_f(event_indx{roi}{trial}, trial);
   end
end

maxT = max(dFF_time(:));
minT = min(dFF_time(:));
histbins = (minT:corr_bin:maxT);
timehist = nan(size(histbins,2),nROI);
parfor roi=1:nROI
    evtm = [];
    for trial = 1:nTrials
        evtm = [evtm; event_times{roi}{trial} ];
    end
    timehist(:,roi) = hist(evtm, histbins);
end
maxlag = 5;
corrs = nan(2*maxlag+1, nROI, nROI);
pkcorr = nan(nROI, nROI);
zerocorr = pkcorr;
for roi1 = 1:nROI
    for roi2 = 1:nROI
%        [corrs(:,roi1, roi2), lags] = xcorr( timehist(:,roi1) - mean(timehist(:,roi1)), timehist(:,roi2) - mean(timehist(:,roi2)), maxlag, 'coeff');
       [corrs(:,roi1, roi2), lags] = xcorr( timehist(:,roi1) , timehist(:,roi2), maxlag, 'coeff');
        pk = max( abs(corrs(:,roi1,roi2)));
       pks = corrs( abs(corrs(:,roi1,roi2))==pk, roi1, roi2 );
       pkcorr(roi1, roi2) = pks(1);
       zerocorr(roi1,roi2) = corrs(maxlag+1, roi1,roi2);
    end
end

%%% between calcium traces
maxlag = 5;
% corrs2 = nan(2*maxlag+1, nROI, nROI
corrs2 = nan(nROI, nROI);
pval2=corrs2;
pkcorr2 = nan(nROI, nROI);
zerocorr2 = pkcorr;
for roi1 = 1:nROI
    for roi2 = 1:nROI
%        [corrs2(:,roi1, roi2), lags2] = xcorr( reshape(temp_dff(:,:,roi1),ceil(Tm/DS_rate)*nTrials,1) - mean(reshape(temp_dff(:,:,roi1),ceil(Tm/DS_rate)*nTrials,1)), reshape(temp_dff(:,:,roi2),ceil(Tm/DS_rate)*nTrials,1) - mean(reshape(temp_dff(:,:,roi2),ceil(Tm/DS_rate)*nTrials,1)), maxlag, 'coeff');
       [corrs2(roi1, roi2), pval2(roi1,roi2)] = corr( reshape(temp_dff(:,:,roi1),ceil(Tm/DS_rate)*nTrials,1) , reshape(temp_dff(:,:,roi2),ceil(Tm/DS_rate)*nTrials,1), 'type', 'Spearman');
%        pk = max( abs(corrs2(:,roi1,roi2)));
%        pks = corrs2( abs(corrs2(:,roi1,roi2))==pk, roi1, roi2 );
%        pkcorr2(roi1, roi2) = pks(1);
%        zerocorr2(roi1,roi2) = corrs2(maxlag+1, roi1,roi2);
    end
end
zerocorr2=corrs2;

figure()
imagesc(pkcorr, [-1,1])
title('Peak correlations for calcium events')

figure()
imagesc(zerocorr, [-1,1])
title('Zero lag correlations for calcium events')

% 
% figure()
% imagesc(pkcorr2, [-1,1])
% title('Peak correlations for calcium traces')

figure()
imagesc(zerocorr2, [-1,1])
title('Zero lag correlations for calcium traces')

end

