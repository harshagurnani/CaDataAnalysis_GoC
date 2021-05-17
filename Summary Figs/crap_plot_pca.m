nROI=p.nROIs;
for sn=1:nSess
    sess = session_types{use_sessions(sn)};
for roi=1:nROI
    [WPOI.(sess).n{roi}, WPOI.(sess).t{roi}] = get_TSOI(Norm.(sess).n(:,roi), Norm.(sess).t(:,roi),Whisk.Detected.(sess),[200 500]);
end
[WPOI.(sess).wMI, WPOI.(sess).wMIt] = get_TSOI( VidMI.(sess).(Whisk.vidseg)(:,1), VidMI.(sess).(Whisk.vidseg)(:,2), Whisk.Detected.(sess), [200 500]);
nWPOI = length(WPOI.(sess).n{roi});
figure; hold on;
allt=0;
npoi_max=30;
for bn=1:nWPOI
    if mod(bn,npoi_max)==1, allt=0; end
   tshift = WPOI.(sess).t{1}{bn}(1)-allt;
   for roi=1:nROI
    plot(WPOI.(sess).t{roi}{bn}-tshift, WPOI.(sess).n{roi}{bn}-prctile(WPOI.(sess).n{roi}{bn},5)+nROI*floor(bn/npoi_max)-0.7*roi, 'k')
   end
   plot(WPOI.(sess).wMIt{bn}-tshift, smooth(WPOI.(sess).wMI{bn},30)*2+nROI*floor(bn/npoi_max), 'g') 
   allt = allt+(WPOI.(sess).t{1}{bn}(end)-WPOI.(sess).t{1}{bn}(1))+1000;
end
title(sess)
end



% plot resp of PCA and whisk

for sn=1:nSess
    sess = session_types{use_sessions(sn)};
    nROI=size(PCA.all.soma.proj,1);
for roi=1:nROI
    [WPOI.(sess).n{roi}, WPOI.(sess).t{roi}] = get_TSOI(Norm.(sess).n(:,roi), Norm.(sess).t(:,roi),Whisk.Detected.(sess),[200 500]);
end
[WPOI.(sess).wMI, WPOI.(sess).wMIt] = get_TSOI( VidMI.(sess).(Whisk.vidseg)(:,1), VidMI.(sess).(Whisk.vidseg)(:,2), Whisk.Detected.(sess), [200 500]);
nWPOI = length(WPOI.(sess).n{roi});
figure; hold on;
allt=0;
npoi_max=30;
for bn=1:nWPOI
    if mod(bn,npoi_max)==1, allt=0; end
   tshift = WPOI.(sess).t{1}{bn}(1)-allt;
   for roi=1:nROI
    plot(WPOI.(sess).t{roi}{bn}-tshift, WPOI.(sess).n{roi}{bn}-prctile(WPOI.(sess).n{roi}{bn},5)+nROI*floor(bn/npoi_max)-0.7*roi, 'k')
   end
   plot(WPOI.(sess).wMIt{bn}-tshift, smooth(WPOI.(sess).wMI{bn},30)*2+nROI*floor(bn/npoi_max), 'g') 
   allt = allt+(WPOI.(sess).t{1}{bn}(end)-WPOI.(sess).t{1}{bn}(1))+1000;
end
title(sess)
end
