% explained variance by first component
pcvar = [];

figure; hold on;
for animal=1:nAnimals
for session = 1:nSess
    try
    pcvar = [pcvar; [SummaryData{animal,session}.PCA.all.soma.explained_var(1), SummaryData.PCA{animal}{session}.baseline.n.explained_var(2)-SummaryData.PCA{animal}{session}.baseline.n.explained_var(1)] ];
    plot([1 2], [SummaryData.PCA{animal}{session}.baseline.n.explained_var(1), SummaryData.PCA{animal}{session}.baseline.n.explained_var(2)-SummaryData.PCA{animal}{session}.baseline.n.explained_var(1)] ,'o','markerfacecolor','w','color','w','linewidth',1)
    catch
    end
    
end
end
bb=bar(1,nanmean(pcvar(:,1)));bb.FaceAlpha=0.5;bb.FaceColor=[0 1 0];
errorbar(1,nanmean(pcvar(:,1)),nanstd(pcvar(:,1)),'g','marker','.','linewidth',3)

bb=bar(2,nanmean(pcvar(:,2)));bb.FaceAlpha=0.5;bb.FaceColor=[0 1 0];
errorbar(2,nanmean(pcvar(:,2)),nanstd(pcvar(:,2)),'g','marker','.','linewidth',3)
