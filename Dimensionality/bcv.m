function [expvar] = bcv( data, frac, nIters, maxK )

    [nT, nX] = size(data);
    nRow = floor(nT*frac(1));
    nCol = floor(nX*frac(2));
    expvar = nan( maxK, nIters );
    for iter = 1:nIters
        t1 = randperm(nT);
        x1 = randperm(nX);
        
        test = data( t1(1:nRow), x1(1:nCol) );
        train = data( t1(nRow+1:end), x1(nCol+1:end) );
        x2 = data( t1(1:nRow), x1(nCol+1:end) );
        x3 = data( t1(nRow+1:end), x1(1:nCol) );
        
        [u,s,v] = svd(train);
        for k = 1:maxK
           x4 = u * s(:, 1:k) * v(:, 1:k)' ;
           err = test - x2 * pinv(x4) * x3;
           expvar(k, iter) = 1 - nanmean(var(err,1)./var(test,1));
        end
        
    end


end