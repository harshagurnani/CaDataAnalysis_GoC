%% Get all correlations
% 
clear correlations pw_distance
for jj=1:nFiles
L = load( all_analysis{jj}, 'allAnalysed');
allpcs{jj} = L.allAnalysed.PCA;
end

SummaryData= struct('pca', allpcs );

%%

eigDist = [];
eigNorm = [];

for jj=1:nFiles
   eigval = SummaryData(jj).pca.dff.all.eig_val; 
   if length(eigval)<20
       vals = [eigval; zeros(20-length(eigval),1)];
   else
       vals = eigval(1:20);
   end
   eigDist = [ eigDist; vals'];
   eigNorm = [ eigNorm; vals'/vals(1) ];
end

%% Plotting
% 
% figure; hold on;
% for jj=1:nFiles
%     plot( eigDist(jj,:), 'color', [.2,.2,.2], 'linewidth', .3);
% end
% plot(nanmean(eigDist,1), 'r', 'linewidth', 1.5);
% 


figure;
subplot(2,1,2); hold on;
for jj=1:nFiles
    plot((eigNorm(jj,:)), 'color', [.2,.2,.2], 'linewidth', .1);
end
plot((nanmean(eigNorm,1)), 'r', 'linewidth', 1.5);
xticks([0,10,20])
yticks([0, .5,1])

subplot(2,1,1); hold on;
for jj=1:nFiles
    plot( log(eigNorm(jj,:)), 'color', [.2,.2,.2], 'linewidth', .1);
end
plot(log(nanmean(eigNorm,1)), 'r', 'linewidth', 1.5);
xticks([0,10,20])
yticks([-5,0])
