% plot resp of PCA and whisk

%%
allPCt = []; maxt=0;
for jj=1:nSess
    sess = session_types{use_sessions(jj)};
    totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
    allPCt = [allPCt; Norm.(sess).t + maxt-totalt];
end

allwhisk = []; maxt=0;
alldetectedwhisk = [];
for jj=1:nSess
    sess = session_types{use_sessions(jj)};
    totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
    temp = VidMI.(sess).(Whisk.vidseg)+[0, maxt-totalt]; 
    temp2 = Whisk.Detected.(sess)+maxt-totalt; 
    allwhisk = [allwhisk; temp];
    alldetectedwhisk = [alldetectedwhisk; temp2];
end


extra = [500 1500];
nROI=size(PCA.all.soma.proj,1);
for roi=1:nROI
    [WPOI.pc.n{roi}, WPOI.pc.t{roi}] = get_TSOI(PCA.all.soma.proj(roi,:), allPCt(:,1),alldetectedwhisk,extra);
end
[WPOI.pc.wMI, WPOI.pc.wMIt] = get_TSOI( allwhisk(:,1), allwhisk(:,2), alldetectedwhisk, extra);
nWPOI = length(WPOI.pc.n{roi});


% 

cm=colormap(jet(nROI));
figure; hold on;
allt=0;
npoi_max=30;
for bn=1:nWPOI
   if mod(bn,npoi_max)==1, allt=0; end
   tshift = WPOI.pc.t{1}{bn}(1)-allt;
   for roi=1:nROI
    sub=((prctile(PCA.all.soma.eigvec(:,roi),20)>-0.05))*prctile(WPOI.pc.n{roi}{bn},10) +     ((prctile(PCA.all.soma.eigvec(:,roi),20)<-0.05))*prctile(WPOI.pc.n{roi}{bn},90);  
    plot(WPOI.pc.t{roi}{bn}-tshift, ((prctile(PCA.all.soma.eigvec(:,roi),20)>-0.05)*2-1)*(WPOI.pc.n{roi}{bn}-sub)*0.2+nROI*floor(bn/npoi_max)-0.7*roi-1, 'color', cm(roi,:))
   end
   plot(WPOI.pc.wMIt{bn}-tshift, smooth(WPOI.pc.wMI{bn},30)*2+nROI*floor(bn/npoi_max), 'g') 
   allt = allt+(WPOI.pc.t{1}{bn}(end)-WPOI.pc.t{1}{bn}(1))+1000;
end
title('PC responses to whisking')
