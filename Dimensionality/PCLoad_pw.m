%Pairwise corr of PC Loadings

nPC = size(tmp.proj, 1 );
nROI = p.nROIs;
nr = floor(sqrt(nPC))+1;
pcload_corr = nan( nROI, nROI, nPC );

for jj=1:nPC 
    pcload_corr(:,:,jj) = tmp.eigvec(:,jj)*tmp.eigvec(:,jj)'*tmp.eig_val(jj);
end


%% plot distance wise relationships
pc_distmap = cell(nPC,1);
% 
figure; hold on;
cmap = colormap(jet(nPC));
for jj=1:nPC
    subplot(nr,nr, jj); hold on;
for roi1 = 1:nROI
    for roi2 = roi1:nROI
        xd = p.FOV * nanmean(ROI.Coordinates.ROIs_X(roi1,:) - ROI.Coordinates.ROIs_X(roi2,:))/512;
        yd = p.FOV * nanmean(ROI.Coordinates.ROIs_Z(roi1,:) - ROI.Coordinates.ROIs_Z(roi2,:))/512;
        zd = nanmean(ROI.Coordinates.ROIs_Z(roi1,:) - ROI.Coordinates.ROIs_Z(roi2,:));
        
        currdist = norm( [xd yd zd]);
%         for jj=1:nPC
        pc_distmap{jj} = [pc_distmap{jj}; currdist, pcload_corr( roi1, roi2, jj)] ; 
        if abs(pcload_corr( roi1, roi2, jj))>0.01
        scatter( currdist, pcload_corr( roi1, roi2, jj), 10, 'filled', 'markerfacecolor', cmap(jj,:))
        else
        scatter( currdist, pcload_corr( roi1, roi2, jj), 5, 'markeredgecolor','k')
        end    
    end 
end
end