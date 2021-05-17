%% Lag correlation with behaviour :
%-----------------------------------------
%%All neurons
[~,zid] = sort(p.POI_coordinates.ROIs_Z(:,1));

% shuffle correlations
box = [prctile(Corr.beh.(sess).(beh).conf_interval(:,1),5),  prctile(Corr.beh.(sess).(beh).conf_interval(:,3),95)]; % range of shuffle correlations

figure;
cjet = colormap( jet(nROI) );

xd=Corr.beh.(sess).(beh).lags;
ff=fill([xd(1) xd(1) xd(end) xd(end)]/1000, [flip(box) box],'k');  % shuffle level
hold on; ff.FaceColor = 'w';ff.FaceAlpha=0.2;

for jj=1:nROI,plot(xd/1000, Corr.beh.(sess).(beh).corr(zid(jj),:),'color',cjet(jj,:),'linewidth',2); end               % corrn of individual roi
legend(ff,'Shuffle level')
title( ['Correlation with ', beh, ' in ', sess, ' session'] )
savefig([exp.analysed, '\corr_with_',beh,'_',sess,'.fig'])

%% Correlation of top somatic PCs
figure; colormap(jet);
plot( PC_Corr.(sess).(beh).lags/1000, PC_Corr.(sess).(beh).corr, 'linewidth',2);
title( ['Correlation of somatic PCs with ', beh, ' in ', sess, ' session'])
savefig([exp.analysed, '\soma_pcs_corr_with_',beh,'_',sess,'.fig'])


%% Correlation of top background PCs

figure; colormap(jet);
nPC = min(10, size( PC_Corr.(sess).background.(beh).corr, 1));
plot( PC_Corr.(sess).(beh).lags/1000, PC_Corr.(sess).background.(beh).corr(1:nPC,:), 'linewidth',2);
title( ['Correlation of background PCs with ', beh, ' in ', sess, ' session'])
savefig([exp.analysed, '\background_pcs_corr_with_',beh,'_',sess,'.fig'])

