FracTest = 0.7;
FracNeuX = 0.8;
% keepN = (nanmax(zscore(activity_all.soma),[],1)>4);
keepN = false(1,p.nROIs); keepN(gid) = true;
nROI = sum(keepN);
keepNroi = find(keepN);
T = randperm(size(activity_all.t,1)); trainT = floor(length(T)*FracTest);
N = randperm(nROI);  trainN =floor(length(N)*FracNeuX);

num_dim = 1:floor(0.7*nROI);
nIters =30;
[res, res_train] = deal(nan(nROI-trainN, length(num_dim), nIters));
for jj=1:nIters
T = randperm(size(activity_all.t,1)); 
N = randperm(nROI); 

Neuron_x = activity_all.soma(:,keepNroi(N(1:trainN)))';
Neuron_y =  activity_all.soma(:,keepNroi(N(trainN+1:end)))';
train_ind = T(1:trainT);
[ num_dim, Test_t, y_est, res(:,:,jj), U_x, U_y, y_est_train, res_train(:,:,jj) ] = ...
            peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );
while nanmean(res(1,:,jj),2)<0.2 
T = randperm(size(activity_all.t,1)); 
N = randperm(nROI); 

Neuron_x = activity_all.soma(:,N(1:trainN))';
Neuron_y =  activity_all.soma(:,N(trainN+1:end))';
train_ind = T(1:trainT);
[ num_dim, Test_t, y_est, res(:,:,jj), U_x, U_y, y_est_train, res_train(:,:,jj) ] = ...
            peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );

end
end
PCA.combined.cval.res = res;
PCA.combined.cval.num_dim = num_dim;
PCA.combined.cval.res_train = num_dim;
% save( [exp.analysed, '\PCA_results.mat'], 'PCA', '-v7.3');
figure; hold on;
for jj=1:nIters
plot(nanmean(res(:,:,jj),1));    
end
%%
plot(nanmean(nanmean(res,1),3),'k','linewidth',2);hold on;
ylim([0 1])