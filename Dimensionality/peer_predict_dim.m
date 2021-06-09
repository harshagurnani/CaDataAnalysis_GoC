function [ num_dim, Test_t, y_est, res, U_x, U_y, y_est_train, res_train ] = peer_predict_dim( Neuron_x, Neuron_y, train_ind, num_dim )
%% Dimensionality estimation by Cross-validation for Peer Prediction
% Inputs:
%
%   X - 2D array : [ROI x time]
%       Activity of training neurons
%   Y - 2D array : [ROI x time]
%       Activity of to-be predicted neurons
%   train_ind - Array of indices
%       Which timepoint indices to use for training (between 1 to length of X)
%       If empty, is taken as first half of X
%   num_dim - Row vector 
%       What dimensional subspace to test quality of prediction for
%
% Outputs:
%   num_dim - Row_vector
%       What dim subspaces were tested
%   Test_t - Array of indices
%       Which timepoint indices in Neuron_y were used as test data
%   y_est - Cell array, same size as num_dim
%       Estimated activity of y in testing set for the corresponding dim of
%       the predictor matrix.
%   Res - 2D arrayL [ROI x nDim] = [#ROI_Y x length(num_dim)]
%       Residuals for each neuronY and num_dim tested = var( y_est - y)
%   Ux - Projections of x to new basis
%   Uy - Projections of y to new basis
%   y_est_train and res_train - same as y_est, except estimate for training
%   half to judge "quality of decomposition"
%
%   Author: Harsha G.
%   08-10-2018

% Modified by Alex 2019/07/01

%   Logic of the code:
%   F_1 = [F_x1; F_y1] = Activity in training half
%   F_2 = [F_x2; F_y2] = Activity in test half
%   Aim: Try to predict F_y2 from F_x2 by performing linear reg on F_1
%       F_1 = U * S_1 *v_1'   where U = [U_x; U_y] -> U learnt on training
%       data                  ---> 1
%       F_x1 = U_x * S_1 *v_1' and F_y1 = U_y * S_1 *v_1'
%       V has the trajector of the 'm' components. If we take only k
%       components in U, S, V, we get a reduced rank approximation.
%
%   For test data:
%       F_x2 = U_x * S_2 *v_2' and F_y2(est) = U_y * S_2 *v_2' --->2
%       Estimate S_2 * V_2' as (U_x * (U_x'*U_x)^-1 )' * F_x2
%       and calculate estimate for F_y2 as in eq2.
%       If we only take k columns of U, we get a k-rank approx.

    [nX, nT] = size( Neuron_x);
    [nY, nT2] = size(Neuron_y);
    if nT ~= nT2, error( 'Neuron_x and Neuron_y must have same number of timepoints'); end
        
    % ugly indexing
    allT = 1: nT; Train_t = false( 1, nT ); Test_t = true( 1, nT );
    if isempty(train_ind)
        disp('No value given for training timepoints. Using first half of the array as training data')
        train_ind = 1:floor(nT/2);
    end
    Train_t( train_ind ) = true;
    Test_t( train_ind ) = false;
    
    x_train = Neuron_x(:, Train_t); x_test = Neuron_x( :, Test_t);
    y_train = Neuron_y(:, Train_t); y_test = Neuron_y( :, Test_t);
    
    % Estimate component directions
    [U, ~, ~] = svd( [x_train; y_train], 'econ' );
    U_x = U(1:nX,:);    U_y = U(nX+1:end,:);
    
    
    % Predict for different subspaces of the decomposition
    if isempty(num_dim)
        num_dim = 1:size(U_x,2); 
    end
    nIters = length(num_dim);
    if iscolumn(num_dim)
        num_dim = num_dim'; 
    end
    [y_est_train, y_est] = deal(cell( nIters, 1)); 
    [res, res_train] = deal(nan( nY, nIters));
    
    jj = 1; stop = 0;
    while jj <= nIters && stop == 0
        ndim = num_dim(jj);
        y_est{jj} = U_y(:, 1:ndim) * ( U_x(:, 1:ndim) / (U_x(:, 1:ndim)'*U_x(:, 1:ndim)) )' * x_test;
        y_est_train{jj}  = U_y(:, 1:ndim) * ( U_x(:, 1:ndim) / (U_x(:, 1:ndim)'*U_x(:, 1:ndim)) )' * x_train;
        res(1:nY, jj) = 1 - (var( y_test-y_est{jj}, [],2))./var(y_test,[],2);
        res_train(1:nY, jj) = 1 - (var( y_train-y_est_train{jj}, [],2))./var(y_train,[],2);
        if mean(res_train(:,jj)) < 0
            stop = 1;
        end
        jj = jj + 1;
    end
    
end