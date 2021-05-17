%% Raw spike correlation


nAnimals =  size(SummaryData,1); %1;%
nSess = size(SummaryData,2); %1;%


SpikeCorrs = cell(3,1);
allDist = cell(3,1);
allCorrSig = cell(3,1);
allCorrCI = cell(3,1);

allCorrRate = [5, 10, 20];

corrPar.maxlag = 200; %ms
corrPar.extratime = 0;
corrPar.minPer = [];
corrPar.nShuffle = 200;
corrPar.plot = false;
corrPar.smoothdt = 200;

success = false( nAnimals, nSess, length(allCorrRate) );

ctr=0;
for an =1:nAnimals
for sess = 1:nSess
if ~isempty( SummaryData{an, sess})
    ctr=ctr+1
    nTm = length(SummaryData{an, sess}.FR);
    spikesT = SummaryData{an, sess}.t(1:nTm,SummaryData{an, sess}.usedROI);
    resFR = subspace_svd( SummaryData{an, sess}.FR', -1 )';
    for rr = 1:length(allCorrRate)
        try
        [lags, corrs, ~, conf_interval, sig ] = ...
            all_period_lagcorr_pairwise_dff( [0 spikesT(end)], resFR, spikesT, ...
                         allCorrRate(rr), corrPar.maxlag, corrPar.extratime, corrPar.minPer, corrPar.nShuffle, corrPar.plot  );
        nc = ceil(size(corrs,3)/2);
        corrs = corrs(:,:,nc);
        allid = ~isnan(corrs);
        corrs = corrs( allid );
        ci = reshape( conf_interval, [size(conf_interval,1)*size(conf_interval,1), 4]);
        ci = ci(allid(:), :);
        sig = ( corrs< ci(:,1)) | ( corrs> ci(:,3));
        SpikeCorrs{rr} = [SpikeCorrs{rr}; corrs];
        allDist{rr} = [allDist{rr}; SummaryData{an,sess}.ROI.Distance.pw(allid)];
        allCorrCI{rr} = [allCorrCI{rr}; ci( allid(sig), :) ];
        allCorrSig{rr} = [allCorrSig{rr}; sig ];
        close all
        success( an, sess, rr) = true;
        catch
         success( an, sess, rr) = false;
        end
    end
end    
end
end

