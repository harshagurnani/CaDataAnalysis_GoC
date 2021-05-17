function[ P ]  = test_pca_on_fr( x, xt )
%% find principal components from population firing rate vector.

    [P.zscored_x, P.x_mean, P.x_sig] = zscore(x);  %each column is one timeseries and is now standardized.

    x_cov = cov(P.zscored_x);
    [P.eigvec, P.eig_val ] = eig(x_cov, 'vector');

    % sort based on higest eigenvalues
    [ P.eig_val, eig_id] = sort( abs(P.eig_val), 'descend' );
    P.eigvec = P.eigvec(:, eig_id');

    nroi = size(x,2);
    explained_var = nan(nroi,1);
    for jj=1:nroi,  explained_var(jj) = sum(P.eig_val(1:jj)); end
    P.explained_var = explained_var/explained_var(end);   %percentage variance explained

%     figure; hold on; plot(1:nroi, 0.9*ones(1,nroi),'k--'); plot(1:nroi, P.explained_var, 'marker','o', 'markerfacecolor', 'k'); ylim([0.05 1.05])

    ndim = find(P.explained_var>0.9,1);
    disp(sprintf('Number of components required to explain 90 percent of the variance is: %d', ndim))

    P.proj = P.eigvec(:,1:ndim)' * P.zscored_x';
    P.proj_orig = P.eigvec(:,1:ndim)' * x';
    P.proj_all = P.eigvec' * P.zscored_x';
    %plot the first few eigenvectors
    for jj = 1:min(ndim,10)
        plot( xt(:,1), P.proj(jj,:)+min(ndim,10)-jj);hold on

    end

end