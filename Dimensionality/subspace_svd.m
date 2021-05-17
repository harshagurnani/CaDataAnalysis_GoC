function [ x_new ] = subspace_svd( x, dim_to_keep )
%% Dimensionality estimation by Cross-validation for Peer Prediction
% Inputs:
%
%   x - 2D array : [ROI x time]
%       Activity of training neurons
%   
%   dim_to_keep - Row vector 
%       What dimensions (subspace) to keep. If instead you want to remove,
%       give it as negative
%
% Outputs:
%   x_new - Cell array, same size as num_dim
%       Estimated activity of y in testing set for the corresponding dim of
%       the predictor matrix.
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

    [nX, nT] = size( x );
   
    % Estimate component directions
    [U, S, V] = svd( x, 'econ' );
   
    all_dim = size(U,2);
     
    if any(dim_to_keep<0), use_dim = true(1,all_dim); use_dim(abs(dim_to_keep))=false;
    else,                  use_dim = false(1,all_dim);use_dim(dim_to_keep) = true;
    end
    
    x_new = U(:, use_dim)*S(use_dim,use_dim)*V(:,use_dim)';
    
end