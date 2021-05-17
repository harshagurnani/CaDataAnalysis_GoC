% explained variance by first component
PC.ExpVar.Same = [];
PC.ExpVar.CV = [];
for animal=1:nAnimals
for session = 1:nSess
    
    try
    if SummaryData{animal,session}.p.nROIs>5%15
    pc1_2=[SummaryData{animal,session}.PCA.all.soma.explained_var(1),SummaryData{animal,session}.PCA.all.soma.explained_var(2)-SummaryData{animal,session}.PCA.all.soma.explained_var(1)];
    PC.ExpVar.Same = [PC.ExpVar.Same;  pc1_2];
%     disp([num2str(animal), ' ', num2str(session), ' ', num2str(pc1_2(1))])
    res =nanmean(nanmean(SummaryData{animal,session}.PCA.combined.cval.res,1),3);
%     disp([num2str(animal), ' ', num2str(session), ' ', num2str(res(1)-pc1_2(1))])
    pc1_2cv = [res(1), res(2)-res(1)];
    if ~(any(pc1_2cv)<0)
    PC.ExpVar.CV = [PC.ExpVar.CV; pc1_2cv];
    else
    PC.ExpVar.CV = [PC.ExpVar.CV; [NaN NaN]];    
    end
    end
    catch
    end
    
end
end
PC.ExpVar.CV(PC.ExpVar.CV<0)=NaN;

figure; hold on;
for jj=1:length(PC.ExpVar.Same)
   plot([1 2], PC.ExpVar.Same(jj,:) ,'o','markerfacecolor','b','color','b','linewidth',1,'markersize',4) 
   plot([1.3 2.3], PC.ExpVar.CV(jj,:) ,'o','markerfacecolor','k','color','k','linewidth',1,'markersize',4) 
  
   plot([1 1.3], [PC.ExpVar.Same(jj,1) PC.ExpVar.CV(jj,1)], 'k','linewidth',0.3)   
   plot([2 2.3], [PC.ExpVar.Same(jj,2) PC.ExpVar.CV(jj,2)], 'k','linewidth',0.3)
end


bb=bar(1,nanmean(PC.ExpVar.Same(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(2,nanmean(PC.ExpVar.Same(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(1.3,nanmean(PC.ExpVar.CV(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';
bb=bar(2.3,nanmean(PC.ExpVar.CV(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';

clear b
b(1) = bar(NaN, NaN);b(1).FaceAlpha=0.5;b(1).FaceColor=[0 0 1];
b(2) = bar(NaN, NaN);b(2).FaceAlpha=0.5;b(2).FaceColor='k';


errorbar(2,nanmean(PC.ExpVar.Same(:,2)),nanstd(PC.ExpVar.Same(:,2)),'b','marker','.','linewidth',3)
errorbar(1,nanmean(PC.ExpVar.Same(:,1)),nanstd(PC.ExpVar.Same(:,1)),'b','marker','.','linewidth',3)
errorbar(2.3,nanmean(PC.ExpVar.CV(:,2)),nanstd(PC.ExpVar.CV(:,2)),'k','marker','.','linewidth',3)
errorbar(1.3,nanmean(PC.ExpVar.CV(:,1)),nanstd(PC.ExpVar.CV(:,1)),'k','marker','.','linewidth',3)

legend( b, 'Same session', 'Cross-validated')
ylabel('Explained variance')
xlabel('PC')
xticks([1 2])


%% PC correlation with behaviour
PC.Corr.whisk = [];
PC.Corr.whisk = [];
PC.Corr.loco  = [];
PC.Corr.loco  = [];

PC.Corr.lagwhisk = [];
PC.Corr.lagwhisk = [];
PC.Corr.lagloco = [];
PC.Corr.lagloco = [];

PC.Corr.topwhisk = [];
PC.Corr.toploco = [];

ctr=1;
for animal=1:nAnimals
for session = 1:nSess
    
    if ~isempty(SummaryData{animal,session})
    
    PC_all.whisk = SummaryData{animal,session}.PC_Corr.combined.beh.whisk.corr;
    PC_all.loco = SummaryData{animal,session}.PC_Corr.combined.beh.loco.corr;
    alllags = SummaryData{animal,session}.PC_Corr.combined.beh.whisk.lags;
    alllags2 = SummaryData{animal,session}.PC_Corr.combined.beh.loco.lags;
    %zero-lag corr with whisk
    nbins = floor(size(PC_all.whisk,2)/2);
    PC.Corr.whisk=[PC.Corr.whisk; abs(PC_all.whisk(1,nbins+1)), abs(PC_all.whisk(2,nbins+1))];
    % top pc
    pkpc = find(abs(PC_all.whisk(:,nbins+1))==max(abs(PC_all.whisk(:,nbins+1))),1);
    if isempty(pkpc), pkpc = 1; end
    PC.Corr.topwhisk = [PC.Corr.topwhisk; PC_all.whisk(pkpc,nbins+1), pkpc];
    % lag of top pc
    [pk,id] = findpeaks([abs(PC_all.whisk(pkpc,:)) 0],'NPeaks',1);
    PC.Corr.lagwhisk = [PC.Corr.lagwhisk; pk alllags(id) ];
    
    %zero-lag corr with loco for running animals
    if max(abs(smooth(SummaryData{animal,session}.Speed(:,2),50)))*1.2>10 && ~isempty(PC_all.loco)%at least 10 cm/s
%         disp([num2str(animal), ' ', num2str(session)])
        nbins = floor(size(PC_all.loco,2)/2);
        PC.Corr.loco=[PC.Corr.loco; abs(PC_all.loco(1,nbins+1)), abs(PC_all.loco(2,nbins+1))];
        
        pkpc = find(abs(PC_all.loco(:,nbins+1))==max(abs(PC_all.loco(:,nbins+1))),1);
        if isempty(pkpc), pkpc = 1; end
        PC.Corr.toploco = [PC.Corr.toploco; abs(PC_all.loco(pkpc,nbins+1)), pkpc];
    
        [pk,id] = findpeaks([abs(PC_all.loco(pkpc,:)) 0],'NPeaks',1);
        PC.Corr.lagloco = [PC.Corr.lagloco; pk alllags2(id) ];

    else
        PC.Corr.loco=[PC.Corr.loco; NaN, NaN];
        PC.Corr.lagloco=[PC.Corr.lagloco; NaN, NaN];
        PC.Corr.toploco = [PC.Corr.toploco; NaN NaN ];
    end
    

    end
    
end
end


%%
figure; hold on;

for jj=1:length(PC.Corr.topwhisk)
    plot([1 3],[PC.Corr.whisk(jj,1), PC.Corr.loco(jj,1)],'ko','markerfacecolor','k','markersize',4)
    plot([1.3 3.3],[PC.Corr.whisk(jj,2), PC.Corr.loco(jj,2)],'ro','markerfacecolor','r','markersize',4)
    plot([0.7 2.7],[PC.Corr.topwhisk(jj,1), PC.Corr.toploco(jj,1)],'bo','markerfacecolor','b','markersize',4)

end
bb = bar(1, nanmean(PC.Corr.whisk(:,1)), 0.3); bb.FaceColor = 'k'; bb.FaceAlpha = 0.5;
bb = bar(1.3, nanmean(PC.Corr.whisk(:,2)), 0.3); bb.FaceColor = 'r'; bb.FaceAlpha = 0.5;
bb = bar(0.7, nanmean(PC.Corr.topwhisk(:,1)), 0.3); bb.FaceColor = 'b'; bb.FaceAlpha = 0.5;

bb = bar(3, nanmean(PC.Corr.loco(:,1)), 0.3); bb.FaceColor = 'k'; bb.FaceAlpha = 0.5;
bb = bar(3.3, nanmean(PC.Corr.loco(:,2)), 0.3); bb.FaceColor = 'r'; bb.FaceAlpha = 0.5;
bb = bar(2.7, nanmean(PC.Corr.toploco(:,1)), 0.3); bb.FaceColor = 'b'; bb.FaceAlpha = 0.5;

clear b
b(1) = bar(NaN, NaN); b(1).FaceColor = 'b'; b(1).FaceAlpha = 0.5;
b(2) = bar(NaN, NaN); b(2).FaceColor = 'k'; b(2).FaceAlpha = 0.5;
b(3) = bar(NaN, NaN); b(3).FaceColor = 'r'; b(3).FaceAlpha = 0.5;

ylim([0 1])
xticks([1 3])
xticklabels({'Whisking','Movement'})
ylabel('Correlation coefficient')
legend(b, 'Best PC', 'PC 1', 'PC 2')