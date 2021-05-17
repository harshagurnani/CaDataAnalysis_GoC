allAnalysed.ID = allData.ID;
path_to_save = 'D:\Work\OneDrive - University College London\pubs and work\Golgi in vivo imaging\Paper\Datasets\Extra\';
%% Pairwise distances

allAnalysed.pw_distance = squareform(pdist( [nanmean(allData.ROI.Coordinates.X,2), nanmean(allData.ROI.Coordinates.Y,2), nanmean(allData.ROI.Coordinates.Z,2)])); 

%% Pairwise correlations

allAnalysed.params.corr.pw = Corr.pw.params;
nlags = floor(size(Corr.pw.combined.all.corr,3)/2)+1;

% full dff
allAnalysed.correlations.all.pw.corr = Corr.pw.combined.all.corr(:,:,nlags);
x = Corr.pw.combined.all.corr_sig(:,:,nlags); x(isnan(x)) = 0; x = x+x';
allAnalysed.correlations.all.pw.sig = x;
allAnalysed.correlations.all.pw.conf_interval = Corr.pw.combined.all.conf_interval;

%without pc1
allAnalysed.correlations.noPC1.pw.corr = Corr.pw.noPC1.all.corr(:,:,nlags);
x = Corr.pw.noPC1.all.corr_sig(:,:,nlags); x(isnan(x)) = 0; x = x+x';
allAnalysed.correlations.noPC1.pw.sig = x;
allAnalysed.correlations.noPC1.pw.conf_interval = Corr.pw.noPC1.all.conf_interval;

%without pc2
% allAnalysed.correlations.noPC2.pw.corr = Corr.pw.noPC2.all.corr(:,:,nlags);
% x = Corr.pw.noPC2.all.corr_sig(:,:,nlags); x(isnan(x)) = 0; x = x+x';
% allAnalysed.correlations.noPC2.pw.sig = x;
% allAnalysed.correlations.noPC2.pw.conf_interval = Corr.pw.noPC2.all.conf_interval;
% 
% %without pc3
% allAnalysed.correlations.noPC3.pw.corr = Corr.pw.noPC3.all.corr(:,:,nlags);
% x = Corr.pw.noPC3.all.corr_sig(:,:,nlags); x(isnan(x)) = 0; x = x+x';
% allAnalysed.correlations.noPC3.pw.sig = x;
% allAnalysed.correlations.noPC3.pw.conf_interval = Corr.pw.noPC3.all.conf_interval;
% 

%% behavioural correlations





%% pca results

%dff
allAnalysed.PCA.dff.all = PCA.all.soma;
try, allAnalysed.PCA.dff.background = PCA.all.background; catch, end
allAnalysed.PCA.dff.crossval = PCA.combined.cval;

%spikes
try
allAnalysed.PCA.fr.all = PCA.combined.spikes;
allAnalysed.PCA.fr.crossval = PCA.FR.cval;
catch
    
end

%% dimensionality

% for dff
res = allAnalysed.PCA.dff.crossval.res; res(res<0)=NaN;
x = nanmean(nanmean(res,3),1);
dim = find(x>max(x)*.99); dim = dim(end);
allAnalysed.Dimensionality.dff.cv = [ dim, x(dim)];

x = find(allAnalysed.PCA.dff.all.explained_var>.8,1);
allAnalysed.Dimensionality.dff.EV80 = [ x, allAnalysed.PCA.dff.all.explained_var(x)];

cc = allAnalysed.PCA.dff.all.eig_val;
allAnalysed.Dimensionality.dff.Spectral = sum(cc)^2/ sum(cc.^2);

% for inferred events
res=allAnalysed.PCA.fr.crossval.res; res(res<0)=NaN;
x = nanmean(nanmean(res,3),1);
dim = find(x>max(x)*.99); dim = dim(end);
allAnalysed.Dimensionality.fr.cv = [ dim, x(dim)];

x = find(allAnalysed.PCA.fr.all.explained_var>.8,1);
allAnalysed.Dimensionality.fr.EV80 = [ x, allAnalysed.PCA.fr.all.explained_var(x)];

cc = allAnalysed.PCA.fr.all.eig_val;
allAnalysed.Dimensionality.fr.Spectral = sum(cc)^2/ sum(cc.^2);

%%
save( [path_to_save, allData.ID, '.mat'], 'allData', 'allEvents', 'allAnalysed', '-v7.3');