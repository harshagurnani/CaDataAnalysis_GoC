nROI=p.nROIs;
background = []; soma = [];figure;plot(1:nROI, 0.9*ones(1,nROI),'k--');  hold on; maxt=0; allt=[];
for sn=1:nSess
    sess = session_types{use_sessions(sn)};
    
    
    for roi=1:nROI, Norm.(sess).background(:,roi) = inpaint_nans(Norm.(sess).background(:,roi),4);  Norm.(sess).n(:,roi) = inpaint_nans(Norm.(sess).n(:,roi),4);end
    f0b = repmat(Norm.baseline.F0, size(Norm.(sess).t,1), 1);
    f0 = repmat(Norm.(sess).F0, size(Norm.(sess).t,1), 1);

    background = [background; Norm.(sess).background];
    soma = [soma; ((Norm.(sess).n.*f0 - f0)-f0b)./f0b];
    allt = [allt;  Norm.(sess).t + maxt];
    maxt = maxt+Norm.(sess).t(end-1);
    PCA.(sess).background = test_pca_on_fr( Norm.(sess).background, Norm.(sess).t ) ;
    PCA.(sess).n = test_pca_on_fr( Norm.(sess).n, Norm.(sess).t ) ;    
    plot(1:nROI, PCA.(sess).background.explained_var, 'w--');
    plot(1:nROI, PCA.(sess).n.explained_var, 'g'); ylim([0.05 1.05])
end
PCA.all.background =  test_pca_on_fr( background, allt ) ;
PCA.all.soma       = test_pca_on_fr( soma, allt );
a(1)=plot(1:nROI, PCA.all.background.explained_var, 'w--', 'linewidth',2); 
a(2)=plot(1:nROI, PCA.all.soma.explained_var, 'g','linewidth',2);  title('Explained variance of PCs');
legend(a,'Background','Golgi soma')

save( [exp.analysed, '\PCA_results.mat'], 'PCA', '-v7.3')
savefig([exp.analysed, '\PC_explainedVar.fig'])


%% Plot traces
figure;sess = session_types{use_sessions(1)}; nComp =size(PCA.(sess).background.proj,1);
plot(Norm.(sess).t(:,1)/1000, PCA.(sess).background.proj+15*nComp-15*(1:nComp)')
ylim([-10 15*nComp+10])
xlim([0 102])
yticks;hold on
title('Background PCs')
savefig([exp.analysed, '\BackgroundPCs.fig'])


figure;nComp =size(PCA.(sess).n.proj,1);
plot(Norm.(sess).t(:,1)/1000, PCA.(sess).n.proj+10*nComp-10*(1:nComp)')
ylim([-10 10*nComp+10])
yticks([])
hold on
title('Golgi population PCs')
savefig([exp.analysed, '\SomaticPCs.fig'])
disp('............................ Saved PCA Results')

%% Get correlations
allloco = []; allrun = []; allmov = []; maxt = 0; allenct = [];
% With running
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
    totalt = Norm.(sess).t(end-1);
%    loco         = -Speed.(sess)(:,2);                           % inverted encoder trace
%    runspeed     = -Speed.(sess)(:,2); runspeed(runspeed<0) = 0; % only forward run speed
%    movement     = abs(Speed.(sess)(:,2));                       % any movement rectified
   loco         = smooth(-Speed.(sess)(:,2),50);                           % inverted encoder trace
   runspeed     = -Speed.(sess)(:,2); runspeed(runspeed<0) = 0;            % only forward run speed
   runspeed     = smooth(runspeed, 50);
   movement     = smooth(abs(Speed.(sess)(:,2)),50);                       % any movement rectified
   allloco = [allloco; loco];  allrun = [allrun; runspeed];
   allmov = [allmov; movement];    allenct = [allenct; Speed.(sess)(:,1)+maxt ];
   maxt = maxt+totalt;
   nComp_n = size(PCA.(sess).n.proj,1);
   [PC_Corr.(sess).loco.lags, PC_Corr.(sess).loco.corr, PC_Corr.(sess).loco.period, PC_Corr.(sess).loco.conf_interval, PC_Corr.(sess).loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], loco , Speed.(sess)(:,1), ...
                                                            PCA.(sess).n.proj', repmat(Norm.(sess).t(:,1),1,nComp_n), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
                                                        
   [PC_Corr.(sess).runspeed.lags, PC_Corr.(sess).runspeed.corr, PC_Corr.(sess).runspeed.period, PC_Corr.(sess).runspeed.conf_interval, PC_Corr.(sess).runspeed.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], runspeed, Speed.(sess)(:,1), ...
                                                            PCA.(sess).n.proj', repmat(Norm.(sess).t(:,1),1,nComp_n), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
   [PC_Corr.(sess).mov.lags, PC_Corr.(sess).mov.corr, PC_Corr.(sess).mov.period, PC_Corr.(sess).mov.conf_interval, PC_Corr.(sess).mov.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], movement, Speed.(sess)(:,1), ...
                                                            PCA.(sess).n.proj', repmat(Norm.(sess).t(:,1),1,nComp_n), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
                                                     
   nComp_bg = size(PCA.(sess).background.proj,1);
   [PC_Corr.(sess).background.loco.lags, PC_Corr.(sess).background.loco.corr, PC_Corr.(sess).background.loco.period, PC_Corr.(sess).background.loco.conf_interval, PC_Corr.(sess).background.loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], loco, Speed.(sess)(:,1), ...
                                                            PCA.(sess).background.proj', repmat(Norm.(sess).t(:,1),1,nComp_bg), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );                                                     
   [PC_Corr.(sess).background.runspeed.lags, PC_Corr.(sess).background.runspeed.corr, PC_Corr.(sess).background.runspeed.period, PC_Corr.(sess).background.runspeed.conf_interval, PC_Corr.(sess).background.runspeed.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], runspeed, Speed.(sess)(:,1), ...
                                                            PCA.(sess).background.proj', repmat(Norm.(sess).t(:,1),1,nComp_bg), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );                                                     
   [PC_Corr.(sess).background.mov.lags, PC_Corr.(sess).background.mov.corr, PC_Corr.(sess).background.mov.period, PC_Corr.(sess).background.mov.conf_interval, PC_Corr.(sess).background.mov.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], movement, Speed.(sess)(:,1), ...
                                                            PCA.(sess).background.proj', repmat(Norm.(sess).t(:,1),1,nComp_bg), ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );                                                     
                                                        
end

sess = 'all';
nComp_n = size(PCA.all.soma.proj,1);
[PC_Corr.(sess).loco.lags, PC_Corr.(sess).loco.corr, PC_Corr.(sess).loco.period, PC_Corr.(sess).loco.conf_interval, PC_Corr.(sess).loco.corr_sig ] =  ...
                    all_period_lagcorr_beh_and_dff_new( [0 maxt], allloco , allenct, ...
                                                        PCA.(sess).soma.proj', repmat(allt(:,1),1,nComp_n), ...
                                                        PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );

[PC_Corr.(sess).runspeed.lags, PC_Corr.(sess).runspeed.corr, PC_Corr.(sess).runspeed.period, PC_Corr.(sess).runspeed.conf_interval, PC_Corr.(sess).runspeed.corr_sig ] =  ...
                    all_period_lagcorr_beh_and_dff_new( [0 maxt], allrun, allenct, ...
                                                        PCA.(sess).soma.proj', repmat(allt(:,1),1,nComp_n), ...
                                                        PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
[PC_Corr.(sess).mov.lags, PC_Corr.(sess).mov.corr, PC_Corr.(sess).mov.period, PC_Corr.(sess).mov.conf_interval, PC_Corr.(sess).mov.corr_sig ] =  ...
                    all_period_lagcorr_beh_and_dff_new( [0 totalT], allmov, allenct, ...
                                                        PCA.(sess).soma.proj', repmat(allt(:,1),1,nComp_n), ...
                                                        PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );



% With Whisking
if exist('VidMI', 'var') && exist('Whisk', 'var')
allwhisk=[]; allwhiskt = []; maxt = 0;
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   totalt = Norm.(sess).t(end-1);
   nComp_n = size(PCA.(sess).n.proj,1);
   tbin = floor(100/nanmean(diff(VidMI.(sess).(Whisk.vidseg)(:,2))));
   [~,ia,~] = unique(VidMI.(sess).(Whisk.vidseg)(:,2));
   uniqueWhisk = smooth(VidMI.(sess).(Whisk.vidseg),tbin);
%    uniqueWhisk([1:tbin*3,end-tbin*3:end])=0;
   uniqueWhisk = uniqueWhisk(ia,1);
   uniqueWhiskTime = VidMI.(sess).(Whisk.vidseg)(ia,2);
   id = find(uniqueWhiskTime>=Norm.(sess).t(end-1),1); if ~isempty(id), uniqueWhiskTime = uniqueWhiskTime(1:id); uniqueWhisk = uniqueWhisk(1:id); end
   allwhisk = [allwhisk; uniqueWhisk];
   allwhiskt = [allwhiskt; uniqueWhiskTime+maxt];
   maxt = maxt+totalt;
   [PC_Corr.(sess).whisk.lags, PC_Corr.(sess).whisk.corr, PC_Corr.(sess).whisk.period, PC_Corr.(sess).whisk.conf_interval, PC_Corr.(sess).whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], uniqueWhisk, uniqueWhiskTime,...
                                                            PCA.(sess).n.proj', repmat(Norm.(sess).t(:,1),1,nComp_n),  ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
   nComp_bg = size(PCA.(sess).background.proj,1);
   [PC_Corr.(sess).background.whisk.lags, PC_Corr.(sess).background.whisk.corr, PC_Corr.(sess).background.whisk.period, PC_Corr.(sess).background.whisk.conf_interval, PC_Corr.(sess).background.whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], uniqueWhisk, uniqueWhiskTime,...
                                                            PCA.(sess).background.proj', repmat(Norm.(sess).t(:,1),1,nComp_bg),  ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
end

sess='all'; nComp_n = size(PCA.all.soma.proj,1);
[PC_Corr.(sess).whisk.lags, PC_Corr.(sess).whisk.corr, PC_Corr.(sess).whisk.period, PC_Corr.(sess).whisk.conf_interval, PC_Corr.(sess).whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], uniqueWhisk, uniqueWhiskTime,...
                                                            PCA.(sess).soma.proj', repmat(allt(:,1),1,nComp_n),  ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
end   
% With pupil size
if exist('Pupil', 'var')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   nComp_n = size(PCA.(sess).n.proj,1);
   [PC_Corr.(sess).pupil.lags, PC_Corr.(sess).pupil.corr, PC_Corr.(sess).pupil.period, PC_Corr.(sess).pupil.conf_interval, PC_Corr.(sess).pupil.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], Pupil.(sess).longAxis, Vidtime.(sess).ecam(:,2),...
                                                            PCA.(sess).n.proj', repmat(Norm.(sess).t(:,1),1,nComp_n),  ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
   nComp_bg = size(PCA.(sess).background.proj,1);
   [PC_Corr.(sess).background.pupil.lags, PC_Corr.(sess).background.pupil.corr, PC_Corr.(sess).background.pupil.period, PC_Corr.(sess).background.pupil.conf_interval, PC_Corr.(sess).background.pupil.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], Pupil.(sess).longAxis,Vidtime.(sess).ecam(:,2),...
                                                            PCA.(sess).background.proj', repmat(Norm.(sess).t(:,1),1,nComp_bg),  ...
                                                            PC_Corr.params.ds_rate, PC_Corr.params.maxlag, [], [], PC_Corr.params.nshuffle, false );
end 
end
save([exp.analysed, '\PC_Corr.mat'], 'PC_Corr','-v7.3');
disp('............................ Saved PC Correlations')
