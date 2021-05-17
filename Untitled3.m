figure;hold on;
for roi=1:8;plot([1 2 3],[BestPCCorr.whisk(roi,1),BestPCCorr.loco(roi,1),BestPCCorr.pupil(roi,1)],'o-','MarkerFaceColor','w','markersize',4,'color','w'); end
bb=bar(1,nanmean(BestPCCorr.whisk(:,1)));bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
errorbar(1,nanmean(BestPCCorr.whisk(:,1)),nanstd(BestPCCorr.whisk(:,1)),'b','marker','.','linewidth',3)
bb=bar(2,nanmean(BestPCCorr.loco(:,1)));bb.FaceAlpha=0.5;bb.FaceColor=[1 0 0];
errorbar(2,nanmean(BestPCCorr.loco(:,1)),nanstd(BestPCCorr.loco(:,1)),'r','marker','.','linewidth',3)
bb=bar(3,nanmean(BestPCCorr.pupil(:,1)));bb.FaceAlpha=0.5;bb.FaceColor=[0 1 0];
errorbar(3,nanmean(BestPCCorr.pupil(:,1)),nanstd(BestPCCorr.pupil(:,1)),'g','marker','.','linewidth',3)
xlim([0.5 3.5])


figure;hold on;
for roi=1:8;plot([1 2 3],[BestPCCorr.whisk(roi),BestPCCorr.loco(roi),BestPCCorr.pupil(roi)],'o-','MarkerFaceColor','w','markersize',4,'color','w'); end
bb=bar(1,nanmean(BestPCCorr.whisk));bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
errorbar(1,nanmean(BestPCCorr.whisk),nanstd(BestPCCorr.whisk),'b','marker','.','linewidth',3)
bb=bar(2,nanmean(BestPCCorr.loco));bb.FaceAlpha=0.5;bb.FaceColor=[1 0 0];
errorbar(2,nanmean(BestPCCorr.loco),nanstd(BestPCCorr.loco),'r','marker','.','linewidth',3)
bb=bar(3,nanmean(BestPCCorr.pupil));bb.FaceAlpha=0.5;bb.FaceColor=[0 1 0];
errorbar(3,nanmean(BestPCCorr.pupil),nanstd(BestPCCorr.pupil),'g','marker','.','linewidth',3)
xlim([0.5 3.5])


bestpccorr

for roi=1:8
    id=~isnan(TopPCCorr.baseline.whisk(roi,:));
    if any(id)
    id = find(abs(TopPCCorr.baseline.whisk(roi,id))==max(abs(TopPCCorr.baseline.whisk(roi,id))));
    BestPCCorr.whisk(roi) = TopPCCorr.baseline.whisk(roi,id);
    else
         BestPCCorr.whisk(roi) = NaN;
    end
    
    id=~isnan(TopPCCorr.baseline.loco(roi,:));
    if any(id)
    id = find(abs(TopPCCorr.baseline.loco(roi,id))==max(abs(TopPCCorr.baseline.loco(roi,id))));
    BestPCCorr.loco(roi) = TopPCCorr.baseline.loco(roi,id);
    else
         BestPCCorr.loco(roi) = NaN;
    end
    
    id=~isnan(TopPCCorr.baseline.pupil(roi,:));
    if any(id)
    id = find(abs(TopPCCorr.baseline.pupil(roi,id))==max(abs(TopPCCorr.baseline.pupil(roi,id))));
    BestPCCorr.pupil(roi) = TopPCCorr.baseline.pupil(roi,id);
    else
         BestPCCorr.pupil(roi) = NaN;
    end
    
end