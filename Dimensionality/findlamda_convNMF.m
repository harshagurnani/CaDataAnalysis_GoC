X= activity_all.spikes'; %neurons  x time

nNeurons = size(X,1);
nBins = ceil(p.acquisition_rate/10); %100 ms sequence

% to redo qqqqq - chunk into more chunks/maybe session wise?
nTrain = ceil(0.6*size(X,2)); %Partition 60% into training
trainX = X(:, 1:nTrain);
testX = X(:, nTrain+1:end);

%% Set up convNMF size
L = nBins; %max length of sequence
K = 4;%ceil(nNeurons/4);   %number of factors


nLambdas = 20;
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 

loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',L,...
         'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end

%% plot costs as a function of lambda
windowSize = 3;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization);
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs);
R = (regularization-minRs)/(maxRs-minRs);
Cs = filtfilt(b,a,cost);
minCs =  prctile(cost,10); maxCs =  prctile(cost,90);
Cs = (Cs -minCs)/(maxCs-minCs);
C = (cost -minCs)/(maxCs-minCs);
clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%%
lambda = 0.005;
nIter = 20;
K = 4;
for iteri = 1:nIter
    [W, H, ~,loadings(iteri,:),power]= seqNMF(trainX,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testX,W,p);
%     W = W(:,is_significant(iteri,:)==1,:); 
%     H = H(is_significant(iteri,:)==1,:); 
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 300; 
    clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end))
    display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
end