function [ev, ev_train, ev_shuff, ev_train_shuff] = pcr_regression_cv( y, yt, x, xt, fs, lambda2, nCV, blockLen )
%% get crossvalidated explained variance for linear fit of y to the combined predictors in x, after sampling at fs
% y is a behaviour , yt are the corresponding timestamps
% x is a timexNeuron predictor matrix, xt are the corresponding timepoints
% fs is sampling rate (in Hz) where all times were given in ms
% lambda2 is an L2 penalty weight

startT = max( min(xt(:), [], 'omitnan'), min(yt(:), [], 'omitnan') );
lastT  = min( max(xt(:), [], 'omitnan'), max(yt(:), [], 'omitnan') );
nX = size(x,2);

if ~isempty(fs), tm = startT:1000/fs:lastT; 
else
    tm = xt; 
    tm( tm<startT) = []; tm( tm>lastT) = []; 
end

% interpolate all to same queried time (fs Hz)
Yq = interp1( yt, y, tm);
if ~isempty(fs)
   Xq = nan( length(tm), nX );
   for jj=1:nX
      Xq(:,jj) = interp1( xt, x(:,jj), tm ); 
   end
else
   Xq = x;
   Xq( xt<startT | xt>lastT,:) = []; 
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


ev = nan( nX, nCV );
ev_train = nan( nX, nCV );

acqRate = 1000/nanmean(diff(tm)); %blockLen = 0.5;
ind = block_shuffle_time(length(tm),acqRate,blockLen); %randperm( length(tm) );



for jj = 1:nCV
    i1 = (jj-1)*nT_test + 1;
    i2 = min( jj*nT_test, length(tm) );
    trainT = true( length(tm),1 );
    trainT( ind(i1:i2) ) = false;
    
    Ytr = Yq( trainT, :);
    Yte = Yq(~trainT, :);
%     n = sum(isnan(Xtr(:))) + sum(isnan(Xte(:))) + sum(isnan(Ytr(:))) + sum(isnan(Yte(:)));
%     fprintf( 'Converting %d NaN to 0 \n', n )
%     Xtr(isnan(Xtr)) = 0;  Xte(isnan(Xte)) = 0; 
%     Ytr(isnan(Ytr)) = 0;  Yte(isnan(Yte)) = 0; 
    for cc = 1:nX
        Xtr = Xq( trainT,[1,2:cc+1]);   Xte = Xq( ~trainT, [1,2:cc+1]); 
        beta = pinv(Xtr'*Xtr + lambda2*eye(cc+1))*Xtr'*Ytr;

        Ypr = Xte * beta ;
        Ypr_train = Xtr * beta ;

        ss = sum ((Yte - mean(Yte)).^2 );
        ress = sum ((Yte - Ypr ).^2 );
        ev( cc, jj) = 1 - ress/ss;

        ss_train = sum ((Ytr - mean(Ytr)).^2 );
        ress_train = sum ((Ytr - Ypr_train ).^2 );
        ev_train( cc, jj) = 1 - ress_train/ss_train;

    end
%     subplot(2,1,1);plot(Yte); hold on; plot(Ypr+30);hold off;
%     subplot(2,1,2); plot(Ytr); hold on; plot(Ypr_train+30); hold off;
   
end


% shuffle control
nShuff = 10;
ev_shuff = nan( nX, nShuff );
ev_train_shuff = nan( nX, nShuff );
for jj=1:nShuff
    
    ind2 = block_shuffle_time(length(tm),acqRate,blockLen);
    Yq_shuff = Yq( ind2, : );   % shuffle up behavior
    
    i1 = (jj-1)*nT_test + 1;
    i2 = min( jj*nT_test, length(tm) );
    trainT = true( length(tm),1 );
    trainT( ind(i1:i2) ) = false;
    
    Ytr = Yq_shuff( trainT, :);
    Yte = Yq_shuff(~trainT, :);
    
    for cc = 1:nX
        Xtr = Xq( trainT,[1,2:cc+1]);   Xte = Xq( ~trainT, [1,2:cc+1]); 
        beta = pinv(Xtr'*Xtr + lambda2*eye(cc+1))*Xtr'*Ytr;

        Ypr = Xte * beta ;
        Ypr_train = Xtr * beta ;

        ss = sum ((Yte - mean(Yte)).^2 );
        ress = sum ((Yte - Ypr ).^2 );
        ev_shuff( cc, jj) = 1 - ress/ss;

        ss_train = sum ((Ytr - mean(Ytr)).^2 );
        ress_train = sum ((Ytr - Ypr_train ).^2 );
        ev_train_shuff( cc, jj) = 1 - ress_train/ss_train;

    end
    
end
