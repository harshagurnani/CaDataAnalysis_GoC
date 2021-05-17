p2 = p; p2.smooth_causal = true; p2.smooth_scale = 5;
for kk=1:nSess
sess = session_types{use_sessions(kk)};
[newNorm.(sess).n, newF0.(sess)] = normalize_and_smooth( reshape(Seg.(sess).r, [p.nTimepoints, p.nTrials, p.nROIs]), reshape(Seg.(sess).t, [p.nTimepoints, p.nTrials, p.nROIs] ), p2);
newNorm.(sess).n=reshape(newNorm.(sess).n, [p.nTimepoints*p.nTrials, p.nROIs]);
end

dt = nanmean(diff(Norm.baseline.t(:,1)));
dt2 = nanmean(diff(VidMI.baseline.whiskpadR_wcam(:,2)));
xwindow = [-400 2000];

indn = floor(xwindow/dt);
indv = floor(xwindow/dt2);

xt = indn(1)*dt:dt:indn(2)*dt;
xt2 = indv(1)*dt2:dt2:indv(2)*dt2;

wtrace = cell(p.nROIs  ,1);
wmitrace = [];
%%
for kk=1:nSess
sess = session_types{use_sessions(kk)};
for jj=1:size(Whisk.Detected.(sess),1)
t1 = find(Norm.(sess).t(:,1)>Whisk.Detected.(sess)(jj,1), 1);
if ~isempty(t1)
for roi=1:p.nROIs
    if length(Norm.(sess).n)>t1+indn(2) && (t1>-indn(1))
    wtrace{roi} = [wtrace{roi}, newNorm.(sess).n(t1+indn(1):t1+indn(2),roi)];
    end
end
t2 = find(VidMI.(sess).whiskpadR_wcam(:,2)>Whisk.Detected.(sess)(jj,1), 1);
if length(VidMI.(sess).whiskpadR_wcam)>t2+indv(2) && t2>-indv(1)
    wmitrace = [wmitrace, VidMI.(sess).whiskpadR_wcam(t2+indv(1):t2+indv(2),1)];
end
end
end
end

wtrace = reshape(wtrace,[1,1,p.nROIs]);
wtrace = cell2mat(wtrace);

%% 
minr = prctile(nanmean(wtrace,3), 5, 2);
maxr = prctile(nanmean(wtrace,3), 95, 2);
stdr = nanstd(nanmean(wtrace,3),[],2);%1/(size(wtrace,2)-1)*

allr = nanmean(nanmean(wtrace, 3),2);
figure;hold on;
a = polyshape([xt, fliplr(xt)], [minr; flipud(maxr)]');
a = polyshape([xt, fliplr(xt)], [allr-stdr; flipud(allr+stdr)]');
plot(a)
plot(xt,allr )
plot([0 0],[0 0.8],'k--')
xlim([-400 2000])
ylim([0 0.8])
