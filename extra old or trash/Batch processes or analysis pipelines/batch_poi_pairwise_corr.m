function [ corrs, lags, pkperiod_corr, pkcorr, used_poi, corrparams] = batch_poi_pairwise_corr( poi, dff, dff_time, extra, dsrates, maxlag, doplot)
%% Pairwise correlations for all ROI in dff - for different kinds of periods (whisking, locomotion etc) and different timebins (downsampled rate dsrates)
%   - poi:       CELL ARRAY Each cell has a 2-column array of start and end
%                times of different periods of interest by one category.
%                Eg. cell 1 can be locomotion periods, cell 2 can be
%                various whisking periods, cell 3 can be paw movement or
%                any behaviour of interest.
%   - dff:       3-D array (Time x Trial x ROI) of dff (or zcored F)
%   - dff_time:  3-D array (Time x Trial x ROI) of acquisition times
%   - extra:     DOUBLE  Time added before/after each period, and subtracted
%                from its negative "gap" periods
%   - dsrates:   Row vector (1xN) of different downsampled rates (for dff)
%                for calculating correlations= 1/(Time bin for correlation)
%                - Has to be divisor of 1000.
%   - maxlag:    DOUBLE Maximum lag (in ms) for computing correlations
%   - doplot:    LOGICAL Plot all correlation graphs?

nPOI = size(poi,2);
Tm = size(dff,1); nTrials = size(dff,2); nROI = size(dff,3);
last_tm = dff_time(end);
corrparams.extratime = extra;
corrparams.dsrates = dsrates;
corrparams.maxlag = maxlag;

for prd = 1:nPOI
    
     lp = poi{prd};
    % Create its negative - gaps between period of interest
    nLP = size(lp,1);
    nNonLP = nLP+1;
    nonLP = nan(nNonLP, 2);
    nonLP(1,:) = [0, lp(1)];
    for kk=2:nLP
    nonLP(kk,:) =[ lp(kk-1,2), lp(kk,1)];
    end

    nonLP(nNonLP,:) = [lp(nLP,2), last_tm];
    
    gap_poi{prd} = nonLP;
end

corrparams.poi =poi;
corrparams.gap_poi =gap_poi;

for jj = 1:size(dsrates,2)
    corrrate = dsrates(jj);
    
    for prd = 1:nPOI
        
        [lags{jj}{prd}{1}, corrs{jj}{prd}{1}, pkperiod_corr{jj}{prd}{1}, pkcorr{jj}{prd}{1}, used_poi{jj}{prd}{1}] =  pairwise_period_lagcorr(poi{prd}, reshape(dff,Tm*nTrials,nROI), reshape(dff_time,Tm*nTrials,nROI),corrrate, maxlag, extra, doplot );
        [lags{jj}{prd}{2}, corrs{jj}{prd}{2}, pkperiod_corr{jj}{prd}{2}, pkcorr{jj}{prd}{2}, used_poi{jj}{prd}{2}] =  pairwise_period_lagcorr(gap_poi{prd}, reshape(dff,Tm*nTrials,nROI), reshape(dff_time,Tm*nTrials,nROI), corrrate, maxlag, -extra, doplot );
    end
    
    

end

end