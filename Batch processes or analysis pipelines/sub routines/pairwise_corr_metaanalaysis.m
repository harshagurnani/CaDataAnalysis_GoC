maxT = nanmax( activity_all.t(:,1));
activity_all.soma_noPC1 = subspace_svd( activity_all.soma', -1 );
useNaN = [];
for roi=1:p.nROIs
useNaN = [useNaN; find(isnan(activity_all.soma_noPC1(roi,:)'))];    
end
useNaN = unique(useNaN);
activity_all.t_noPC1 = activity_all.t;
activity_all.soma_noPC1(:,useNaN) =[];
activity_all.t_noPC1(useNaN,:) = [];

% normal cross corr
[Corr.combined.pw.all.lags, Corr.combined.pw.all.corr, Corr.combined.pw.all.period, Corr.combined.pw.all.conf_interval, Corr.combined.pw.all.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxT], activity_all.soma, activity_all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, 0, [], Corr.pw.params.nshuffle, false );
                            
% noPC1 cross corr
[Corr.combined.pw.noPC1.lags, Corr.combined.pw.noPC1.corr, Corr.combined.pw.noPC1.period, Corr.combined.pw.noPC1.conf_interval, Corr.combined.pw.noPC1.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxT], activity_all.soma_noPC1', activity_all.t_noPC1, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, 0, [], Corr.pw.params.nshuffle, false );


nPairs = 0.5*(p.nROIs)*(p.nROIs-1);
[CorrAnalysis.allD.rPW, CorrAnalysis.allD.rAP, CorrAnalysis.allD.rML ] = deal( nan(nPairs,1) );
[CorrAnalysis.peakCC, CorrAnalysis.peakCC_noPC1,  CorrAnalysis.lagCC, CorrAnalysis.lagCC_noPC1] = deal( nan(nPairs,1) );
[CorrAnalysis.peakW, CorrAnalysis.peakW_noPC1, CorrAnalysis.peakProm, CorrAnalysis.peakProm_noPC1] = deal( nan(nPairs,1) );
[CorrAnalysis.peakSig, CorrAnalysis.peakSig_noPC1] =  deal( nan(nPairs,1) );


nLags = length(Corr.combined.pw.all.lags);
[CorrAnalysis.allC, CorrAnalysis.allC_noPC1 ] = deal(nan( nLags, nPairs ));

CorrAnalysis.Locations = ROI.Coordinates;
CorrAnalysis.px2um = 460/512;


ctr=0;
for roi1=1:p.nROIs
for roi2=roi1+1:p.nROIs
    ctr=ctr+1;
    
    % CCgram
    CorrAnalysis.allC( :, ctr ) = Corr.combined.pw.all.corr(roi1, roi2,:);
    CorrAnalysis.allC_noPC1( :, ctr ) = Corr.combined.pw.noPC1.corr(roi1, roi2,:);
    
    %Peak CC
    [pk1,pk2,pk3,pk4]= findpeaks(CorrAnalysis.allC( :, ctr ),'NPeaks',1);
    if ~isempty(pk1), CorrAnalysis.peakCC(ctr)= pk1;    CorrAnalysis.lagCC(ctr) = Corr.combined.pw.all.lags(pk2); 
                      CorrAnalysis.peakW(ctr) = pk3*50; CorrAnalysis.peakProm(ctr) = pk4;
                      CorrAnalysis.peakSig(ctr) = Corr.combined.pw.all.corr_sig( roi1, roi2, pk2 );
    end
    
    [pk1,pk2,pk3,pk4]= findpeaks(CorrAnalysis.allC_noPC1( :, ctr ),'NPeaks',1);
    if ~isempty(pk1), CorrAnalysis.peakCC_noPC1(ctr)= pk1;    CorrAnalysis.lagCC_noPC1(ctr) = Corr.combined.pw.noPC1.lags(pk2); 
                      CorrAnalysis.peakW_noPC1(ctr) = pk3*50; CorrAnalysis.peakProm_noPC1(ctr) = pk4;
                      CorrAnalysis.peakSig_noPC1(ctr) = Corr.combined.pw.noPC1.corr_sig( roi1, roi2, pk2 );
    end
    
    %Distances
    CorrAnalysis.allD.rPW(ctr) = ROI.Distance.pw(roi1, roi2);
    CorrAnalysis.allD.rAP(ctr) = abs(ROI.Coordinates.ROIs_X(roi1,1)-ROI.Coordinates.ROIs_X(roi2,1));
    CorrAnalysis.allD.rML(ctr) = abs(ROI.Coordinates.ROIs_Y(roi1,1)-ROI.Coordinates.ROIs_Y(roi2,1));
end
end