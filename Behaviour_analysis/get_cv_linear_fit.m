function [ev, ev_train] = get_cv_linear_fit( y, yt, x, fs,  lambda2, nCV )
%% get crossvalidated explained variance for linear fit of y to the combined predictors in x, after sampling at fs
% y is a timexNeurons , yt are the corresponding timestamps
% x is a cell array of 2-column vectors, x{jj} = [ time_jj, value_jj]_t
% fs is sampling rate (in Hz) where all times were given in ms
% lambda2 is an L2 penalty weight

startT = min( yt(:), [], 'omitnan' );
lastT  = max( yt(:), [], 'omitnan' );
nX = length(x); nNeu = size(y,2);
for jj=1:nX, startT = [lastT, min(x{jj}(:,1), [], 'omitnan')];  lastT = [startT, max(x{jj}(:,1), [], 'omitnan')];  end

startT = min(startT);
lastT = max( lastT);
tm = startT:1000/fs:lastT;

% interpolate all to same queried time (fs Hz)
Yq = nan( length(tm), nNeu );
for jj=1:nNeu
    Yq(:,jj) = interp1(yt(:,jj), y(:,jj), tm );
end
Xq = nan( length(tm), nX );
for jj=1:nX
    Xq(:,jj) = interp1( x{jj}(:,1), x{jj}(:,2), tm );
end
Xq = [ones(length(tm),1), Xq];  %add col of 1

% Remove NaN values
ind1 = find( any( isnan(Xq), 2 ) );
ind2 = find( any( isnan(Yq), 2 ) );
ind = unique( [ind1; ind2] );
Xq(ind,:) = [];
Yq(ind,:) = [];
tm(ind) = [];

% nCV = 10; 
nT_test = floor( length(tm)/nCV );

ev = nan( nNeu, nCV );
ev_train = nan( nNeu, nCV );

ind = randperm( length(tm) );
for jj = 1:nCV
    i1 = (jj-1)*nT_test + 1;
    i2 = min( jj*nT_test, length(tm) );
    trainT = true( length(tm),1 );
    trainT( ind(i1:i2) ) = false;
    
    Xtr = Xq( trainT,:);   Ytr = Yq( trainT, :);
    Xte = Xq( ~trainT, :); Yte = Yq(~trainT, :);
%     n = sum(isnan(Xtr(:))) + sum(isnan(Xte(:))) + sum(isnan(Ytr(:))) + sum(isnan(Yte(:)));
%     fprintf( 'Converting %d NaN to 0 \n', n )
%     Xtr(isnan(Xtr)) = 0;  Xte(isnan(Xte)) = 0; 
%     Ytr(isnan(Ytr)) = 0;  Yte(isnan(Yte)) = 0; 
    
    beta = pinv(Xtr'*Xtr + lambda2*eye(nX+1))*Xtr'*Ytr;
    
    Ypr = Xte * beta ;
    Ypr_train = Xtr * beta ;
    
    for roi = 1:nNeu
       ss = sum ((Yte(:,roi) - mean(Yte(:,roi))).^2 );
       ress = sum ((Yte(:,roi) - Ypr(:,roi) ).^2 );
       ev( roi, jj) = 1 - ress/ss;
       
       ss_train = sum ((Ytr(:,roi) - mean(Ytr(:,roi))).^2 );
       ress_train = sum ((Ytr(:,roi) - Ypr_train(:,roi) ).^2 );
       ev_train( roi, jj) = 1 - ress_train/ss_train;
       
    end
    
    
end