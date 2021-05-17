%% pw corr
APResp = [];
APtraj = [];
APtm =[-120,1200];
figure; hold on; ctr=0;
cmap = colormap(redblue(18));
latency = [];
on_off = [];

for animal=1:nAnimals
for session = 1:nSess
    sess='combined';
    if ~isempty(SummaryData{animal,session})
    try
       allAP = SummaryData{animal, session}.noPC1_resp.AP;
       allN  = squeeze(nanmean(allAP.AlignedN,2));
       allT = allAP.AlignedT;
       nroi = size(allN,2);
       allPC = SummaryData{animal, session}.PCA.all.soma.eigvec;
       for roi=1:nroi
           t2000= find(allT>1400,1);
           t50  = find(allT>-50,1);
           apr = zscore(allN(:,roi));
           [y,b] = resample_ds( apr, allT, 25, APtm);
           APtraj = [APtraj; y];
           apr =apr(t50:t2000);
           tmp = allT(t50:t2000);
           maxt = find( apr ==  max(apr),1);
           mint = find( apr ==  min(apr),1);
           maxresp = apr(maxt); minresp = apr(mint);
           if maxresp>-minresp, APResp = [APResp;   maxresp]; latency = [latency; tmp(maxt)+60]; on_off=[on_off;1];
           else,APResp = [APResp;   minresp]; latency = [latency; tmp(mint)+60]; on_off=[on_off;0];
           end
           if abs(min(allN(1:end-15,roi)))< max(allN(1:end-15,roi)),
           cm = floor(10*max(allN(1:end-15,roi)))+10;
       else
           cm = 9- floor(10*abs(min(allN(1:end-15,roi))));
       end
           plot(allT, allN(:,roi),'color', cmap(cm,:));
       end
       ctr=ctr+nroi;
    catch
        
    end
    
    end
end
end

a=area([0 100],[.8,.8],-.8)
a.FaceColor='m';a.FaceAlpha=0.1;a.EdgeAlpha=0;
% plot([100 100],[-.8 .8],'k')

xlim([-300 1500])
yticks([-.8:.4:.8])