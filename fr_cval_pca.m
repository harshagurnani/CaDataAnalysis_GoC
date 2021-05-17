FracTest = 0.4;
FracNeuX = 0.9;
nROI = length(gid);
T = randperm(size(spikesT,1)); trainT = floor(length(T)*FracTest);
N = randperm(nROI);  trainN =floor(length(N)*FracNeuX);

num_dim = 1:floor(0.7*nROI);
nIters =50;
[res, res_train] = deal(nan(nROI-trainN, length(num_dim), nIters));
min6 = true;
for jj=1:nIters
T = randperm(size(spikesT,1)); 
N = randperm(nROI); 

Neuron_x = FR(:,N(1:trainN))';
Neuron_y =  FR(:,N(trainN+1:end))';
train_ind = T(1:trainT);
[ num_dim, Test_t, y_est, res(:,:,jj), U_x, U_y, y_est_train, res_train(:,:,jj) ] = ...
            peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );
while nanmean(nanmean(res(:,:,jj),2),1)<0.2 
T = randperm(size(spikesT,1)); 
N = randperm(nROI); 

Neuron_x = FR(:,N(1:trainN))';
Neuron_y =  FR(:,N(trainN+1:end))';
train_ind = T(1:trainT);
[ num_dim, Test_t, y_est, res(:,:,jj), U_x, U_y, y_est_train, res_train(:,:,jj) ] = ...
            peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );

end
end
PCA.FR.cval.res = res;
PCA.FR.cval.num_dim = num_dim;
PCA.FR.cval.res_train = num_dim;
% save( [exp.analysed, '\PCA_results.mat'], 'PCA', '-v7.3');
figure; hold on;
for jj=1:nIters
plot(nanmean(res(:,:,jj),1));    
end
plot(nanmean(nanmean(res,1),3),'k','linewidth',2);hold on;
ylim([0 1])