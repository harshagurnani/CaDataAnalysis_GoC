%% 
xx=WhiskData.baseline.alltheta_lp(:,2);
[nc, edges, bins]= histcounts( xx, -360:5:360);
mids = edges(2:end)/2+edges(1:end-1)/2
newm = interp1( vtimes, xx(ia), Norm.baseline.t(:,1));
for roi=1:p.nROIs,
for xval = 1:length(mids)
set_resp(xval, roi) = nanmean(Norm.baseline.n(newm>=edges(xval)&newm<edges(xval+1),roi));
end
end
figure;scatter(newm, Norm.baseline.n(:,30),2,'k','filled')
hold on;
plot(mids, set_resp(:,30))
x_new = subspace_svd( Norm.baseline.n', -1 );
[WTrace.(sess).lags, WTrace.(sess).corrs, ~, WTrace.(sess).shuffLevel, WTrace.(sess).corrsig] =      all_period_lagcorr_beh_and_dff( [0 Norm.(sess).t(end)], vid, vtimes, x_new', Norm.(sess).t, Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [],[],Corr.beh.params.nshuffle );
whisk = 2; yy = alltheta_fast.(sess)(:,whisk)-alltheta_lp.(sess)(:,whisk);
whisk = 2; yy = alltheta_lp.(sess)(:,whisk)-alltheta_slow.(sess)(:,whisk);
[WPhase.(sess).lags, WPhase.(sess).corrs, ~, WPhase.(sess).shuffLevel, WPhase.(sess).corrsig] =      all_period_lagcorr_beh_and_dff( [0 Norm.(sess).t(end)], yy(ia), vtimes, x_new', Norm.(sess).t, Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [],[],Corr.beh.params.nshuffle );


%% plotting
[~,id]=sort(sum(WPhase.baseline.corrs(:,20:22),2));
cdata2 = floor(fliplr(WPhase.baseline.corrs(id,:).*abs(WPhase.baseline.corrsig(id,:)))*250+100);
rc = cmap(cdata2,1);
gc = cmap(cdata2,2);
bc = cmap(cdata2,3);
rc = reshape(rc,[40,41]);
gc = reshape(gc,[40,41]);
bc = reshape(bc,[40,41]);
figure;imagesc(cat(3,rc,gc,bc))
xticks([1:10:41])
xticklabels(num2str(WTrace.baseline.lags(1:10:41)'))
xlim([6,36])
[~,id2]=sort(sum(WPhase.baseline.corrs(:,20:22),2),'descend')
figure;
plot(Norm.baseline.t(:,1), x_new(id2,:)'+[1:p.nROIs]*.4)
plot(Norm.baseline.t(:,1)/1000, x_new(id2,:)'+[1:p.nROIs]*.4)
xlim([-5 205])
hold on;
plot(vtimes/1000, yy(ia)/360*5-1,'k','linewidth',1)