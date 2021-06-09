function [cvcurve, peakdim, ExpVar, ExpVar_train, num_dim] = get_cvpca_curves( data, nIters, FracTrainNeu, FracTrainT, maxDim )
%% Function to get cross-validated explained variance curves with multiple iterations
% Inputs:
%   - data          [Time x Neurons 2D array] 
%           Timeseries of dFF or spiking of all neurons
%   - nIters        [INT]
%           Number of iterations for cv-pca
%   - FracTrainNeu  [FLOAT between 0 and 1]
%           Fraction of all neurons in training set
%   - FracTrainT    [FLOAT between 0 and 1]
%           Fraction of all timepoints in training set
%   - maxDim        [NUMERIC]
%           Either dimensionality to test as integer, or as fraction of
%           total number of ROI
%
%
% Outputs:
%   - ExpVar           [ #test_neurons  x  #pcs x #iterations]
%           cross-validated explained variance curves for multiple test
%           neuron (different per iteration) for all iterations and 
%           #components tested
%   - ExpVar           [ #test_neurons  x  #pcs x #iterations]
%           Explained variance on training timepoints
%   - num_dim       [int array]
%           # pcs tested
%
%   Harsha Gurnani. June 2021
    
    [nT, nROI] = size(data);
    trainT = floor(nT*FracTrainT);
    trainN =floor(nROI*FracTrainNeu);

    if maxDim<=1, maxDim = floor(maxDim*nROI); end
    if maxDim>nROI, maxDim = nROI; end
    
    num_dim = 1:maxDim;
    
    [ExpVar, ExpVar_train] = deal(nan(nROI-trainN, length(num_dim), nIters));
    for jj=1:nIters
        T = randperm(nT); 
        N = randperm(nROI); 

        Neuron_x = data(:,N(1:trainN))';            % training neurons
        Neuron_y =  data(:,N(trainN+1:end))';       % test neurons
        train_ind = T(1:trainT);                    % training timepoints
        [ num_dim, Test_t, y_est, ExpVar(:,:,jj), U_x, U_y, y_est_train, ExpVar_train(:,:,jj) ] = ...
                            peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );
%         while nanmean(ExpVar(1,:,jj),2)<0.05 
%             % poor partitioning -  uncomment if needed
%             T = randperm(nT); 
%             N = randperm(nROI); 
% 
%             Neuron_x = data(:,N(1:trainN))';            % training neurons
%             Neuron_y =  data(:,N(trainN+1:end))';       % test neurons
%             train_ind = T(1:trainT);                    % training timepoints
%             [ num_dim, Test_t, y_est, ExpVar(:,:,jj), U_x, U_y, y_est_train, ExpVar_train(:,:,jj) ] = ...
%                                 peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim );
% 
%         end
    end
    ExpVar(ExpVar<0)=0;
    ExpVar_train(ExpVar_train<0)=0;
    
    % curve is averaged across neurons and iterations
    cvcurve = nanmean(nanmean(ExpVar,1),3);
    maxCV = max(cvcurve);
    peakdim = find( cvcurve==maxCV);
    
    figure; hold on;
    for jj=1:nIters
    plot(nanmean(ExpVar(:,:,jj),1));    
    end
    
    plot(cvcurve,'k','linewidth',2);hold on;
    ylim([0 1])
    
    


end