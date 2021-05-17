function [ res, error, maxK ] = bi_cross_val( X, keep_neu, keep_tm, nIters )
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


    [nNeu, nT ]= size(X);
    maxK = nNeu - keep_neu;
    res = nan( keep_neu, maxK, nIters );
    error = nan( maxK, nIters );
    for jj = 1:nIters
        T1 = false( 1, nT );   
        N1 = false( nNeu, 1);  
        
        test_t = randperm( nT, keep_tm);
        test_n = randperm( nNeu, keep_neu );
        
        T1( test_t ) = true; T2 = ~T1;
        N1( test_n ) = true; N2 = ~N1;
        
        [U2, S, V2 ] = svd( X(N2, T2), 'econ' );
        A = X(N1, T1); B = X(N1, T2); C = X(N2,T1);
        
        for dim = 1:maxK
           D_hat = U2(:,1:dim) * S(1:dim, 1:dim) * V2(:,1:dim )';
           A_hat = B*pinv(D_hat)*C;
           res(:, dim, jj)= 1 - (var( A -  A_hat, [],2))./var(A,[],2);
           error( dim, jj) = norm( A-A_hat);
        end
        
    end
    
    
    
end