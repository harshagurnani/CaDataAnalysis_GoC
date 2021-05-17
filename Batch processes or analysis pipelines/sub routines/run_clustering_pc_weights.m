%% clustering based on weights
nComp =size(PCA.all.soma.proj,1);
allwts = tmp.eigvec(:,1:min(nComp,15));%PCA.all.soma.eigvec(:,1:min(nComp,15));
resp  = arrayfun(@(n) norm(allwts(n,:)),1:p.nROIs);
testNC = 8;
for nC=1:testNC
    [cluster(nC).id, cluster(nC).centroid, cluster(nC).sumD] = kmeans(allwts,nC);
    clustSum(nC) = sum( cluster(nC).sumD );
end
nC = find( clustSum == nanmin(clustSum), 1 );
use_pc = [2 3 4];

%% plot clustering results
figure;cmap = colormap(jet(nC)); 
subplot(2,1,1);hold on;
for jj=1:nC
    neu = (cluster(nC).id == jj);
scatter3( allwts(neu,use_pc(1)), allwts(neu,use_pc(2)), allwts(neu,use_pc(3)), 'filled','markerfacecolor',cmap(jj,:))
end
xlabel(['PC', num2str(use_pc(1))])
ylabel(['PC',num2str(use_pc(2))])
zlabel(['PC', num2str(use_pc(3))])


subplot(2,1,2); hold on;
for jj=1:nC
    neu = (cluster(nC).id == jj);
    scatter3( nanmean(ROI.Coordinates.ROIs_X(neu,:),2)*400/512, nanmean(ROI.Coordinates.ROIs_Y(neu,:),2)*400/512,...
        nanmean(ROI.Coordinates.ROIs_Z(neu,:),2), 'filled','markerfacecolor',cmap(jj,:))
end
xlabel('X')
ylabel('Y')
zlabel('Z')

% savefig([exp.analysed,'\Cluster_ID.fig'])


nCROI = arrayfun(@(c) sum(cluster(nC).id==c),1:nC);
maxCSize = max(10,nanmax(nCROI));

figure;
hold on
% mean patches
masks=Seg.masks; nROI = numel(masks);
colormap(gray);
for jj=1:nC
neu = find(cluster(nC).id == jj);
for kk=1:length(neu)
    nn = neu(kk);
%     subaxis(2*nC,maxCSize,(jj-1)*maxCSize + kk,'SpacingVert',0,'SpacingHoriz',0); 
    imagesc(Seg.mean_image{nn}); yticks([]); xticks([]); 
    set(gca,'XColor',cmap(jj,:),'YColor',cmap(jj,:))
end
end

subaxis(2,1,2,'SpacingHoriz',0); hold on;
for jj=1:nC
    neu = (cluster(nC).id == jj);
    
    plot(activity_all.t(:,1)/1000, nanmean(activity_all.soma(:,neu),2)+(nC-jj)+1, 'color',cmap(jj,:))
end
% beh_all.loco(beh_all.loco(:,2)<-.05,2) = NaN; beh_all.loco(:,2)=inpaint_nans(beh_all.loco(:,2));
beh(1)=plot(beh_all.loco(:,1)/1000, smooth(beh_all.loco(:,2),50)*.2+nC+1,'k','linewidth',1.2)
beh(2)=plot(beh_all.whisk(:,2)/1000, beh_all.whisk(:,1)+2+nC,'b','linewidth',1.2)
yticks([])
% savefig([exp.analysed,'\Cluster_traces.fig'])

% save([exp.analysed,'\clustering.mat'],'cluster','clustSum','-v7.3')