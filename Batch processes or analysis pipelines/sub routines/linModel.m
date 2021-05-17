% parameters
dsrate    = 20;
w_smooth  = 60;
s_smooth  = 200;
puffwin   = [-200 1200];
lag_shift = false;
nROI = p.nROIs;

alpha = 0.1;
nCV   = 5;

maxT.AP = Norm.AP.t(end,1);
maxT.baseline = Norm.baseline.t(end,1);

maxF = arrayfun( @(roi) prctile( Norm.baseline.n(:,roi), 95 ), 1:nROI);
% Normalise behaviour to same units

allw = [ VidMI.AP.whiskpadR_wcam(:,1); VidMI.baseline.whiskpadR_wcam(:,1) ] ;
allw( isnan(allw) ) = prctile(allw,5);
minw = prctile(allw, 1);
maxw = prctile(allw-minw, 99);

alls = -[ Speed.AP(:,2); Speed.baseline(:,2) ];
% alls(alls <-0.1)  = NaN; alls = inpaint_nans( alls );
% alls(alls <-0.05) = NaN; alls = inpaint_nans( alls );
% alls(alls <-0.01) = NaN; alls = inpaint_nans( alls );
maxs = prctile(alls, 99);

sess = {'AP', 'baseline'};
for jj=1:2
    sn = sess{jj};    
    [~,idx]= unique(VidMI.(sn).whiskpadR_wcam(:,2));idx = sort(idx,'ascend');
    VidMI.(sn).whiskpadR_wcam=VidMI.(sn).whiskpadR_wcam(idx,:);
    
    [ds_whisk.(sn), ds_t.(sn)] = resample_ds( smooth( VidMI.(sn).whiskpadR_wcam(:,1), w_smooth), VidMI.(sn).whiskpadR_wcam(:,2), dsrate, [0 maxT.(sn)] );
    ds_whisk.(sn) = (ds_whisk.(sn) - minw)/maxw;
    
    spd = -smooth( Speed.(sn)(:,2), s_smooth, 'lowess' ); 
    [tmp, ds_t.(sn)] = resample_ds( spd, Speed.(sn)(:,1), dsrate, [0 maxT.(sn)] );
    tmp(tmp<-.1) = NaN; tmp = inpaint_nans(tmp);
    tmp(tmp<-.05) = NaN; tmp = inpaint_nans(tmp);
    tmp(tmp<-.01) = NaN; tmp = inpaint_nans(tmp);
    ds_locom.(sn) = tmp;%/maxs;
    
end


% Make Puff Matrix
pstart = floor(-puffwin(1)/dsrate);
pend   = floor(puffwin(2)/dsrate);
plen = pend+pstart;
nPuff = size(Puff.AP,1);

ds_puff.baseline = zeros( length(ds_t.baseline), plen );
ds_puff.AP       = zeros( length(ds_t.AP), plen );
for jj=1:nPuff
    t1 = find(ds_t.AP > Puff.AP(jj,1) , 1 );

    for kk=1:plen
    idx = t1-pstart-1+kk;    
    ds_puff.AP(idx,kk) = 1;
    end
    
end

%Normalise DFF
maxF = nan(p.nROIs,1);
for roi=1:p.nROIs
   tmp = arrayfun( @(sn) prctile( Norm.(sess{sn}).n(:,roi), 99.9 ), 1:2 );
   maxF(roi) = max(tmp);
end


%% Fit Sessions Individually
clear B FitInfo
for s1 = 1:2
    sn = sess{s1};
    B.(sn) = cell(nROI,1);
    FitInfo.(sn) = B.(sn);
    
    XNew = [ds_whisk.(sn)', ds_locom.(sn)', ds_puff.(sn) ] ;
    
    for roi=1:nROI
       y = resample_ds( Norm.(sn).n(:, roi) , Norm.(sn).t(:, roi), dsrate, [0 maxT.(sn)]);
       y = 1/maxF(roi) * (y*Norm.(sn).F0(roi) +Norm.(sn).F0(roi) - Norm.baseline.F0(roi))/Norm.baseline.F0(roi);
       
       [B.(sn){roi}, FitInfo.(sn){roi}] = lasso( XNew, y, 'Alpha', alpha, 'CV', nCV); %, 'CV', nCV
    end
    
    
end


%% Cross validate on Baseline from AP
sn = 'baseline';
y = resample_ds( Norm.baseline.n(:, roi) , Norm.baseline.t(:, roi), dsrate, [0 maxT.baseline]);
OldY = nan( length(y), nROI ); NewY = OldY;
XNew = [ds_whisk.(sn)', ds_locom.(sn)', ds_puff.(sn) ] ;

for roi = 1:nROI
       y = resample_ds( Norm.(sn).n(:, roi) , Norm.(sn).t(:, roi), dsrate, [0 maxT.(sn)]);
       OldY(:,roi) = 1/maxF(roi) * (y*Norm.(sn).F0(roi) +Norm.(sn).F0(roi) - Norm.baseline.F0(roi))/Norm.baseline.F0(roi);
       bestidx = FitInfo.AP{roi}.Index1SE;
       behwts = B.AP{roi}(:,bestidx);
       bestidx = FitInfo.(sn){roi}.Index1SE;
       intercept = FitInfo.(sn){roi}.Intercept(bestidx);
       NewY(:,roi) = XNew * behwts +  intercept;
end


%% Plotting
presp = nan(size(B.AP{1},1)-2,nROI);
Wts.AP = nan(nROI,2);
Wts.baseline = Wts.AP;

for s1=1:2
    sn=sess{s1};
    for roi=1:nROI
        bestidx = FitInfo.(sn){roi}.Index1SE;
        Wts.(sn)(roi,:) = B.(sn){roi}(1:2,bestidx);
    end
end

figure;
for roi=1:nROI
   bestidx = FitInfo.AP{roi}.Index1SE; 
   presp(:,roi) = B.AP{roi}(3:end,bestidx);
   
   plot( B.AP{roi}(3:end,bestidx)+roi*.3); hold on;
end

figure;imagesc(presp',[-.3 .3])

