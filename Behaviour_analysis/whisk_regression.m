%% Create predictor matrix

ds_rate = 50;   %Hz
error_frames = [1:84, 2197:2255];

vars_for_ica = {'alltheta_lp', 'alltheta_fast', 'WPos_fast'};
XTmp = [];
n_ica_comp = 3;

for jj=1:length(vars_for_ica)
    XTmp = [XTmp, WhiskData.(sess).(vars_for_ica{jj})];
end
mdl = rica( XTmp, n_ica_comp);
XTmp = transform(mdl, XTmp);

other_vars = {'wpad_slow', 'VidMI_fast'};%
for jj=1:length(other_vars)
   XTmp = [XTmp,  WhiskData.(sess).(other_vars{jj})];
end
XTmp = zscore(XTmp)/5;

XConv = [];
nxconv = 9*size(XTmp,2);
for jj=1:size(XTmp,2)
    for kk=1:9
       XConv = [XConv, conv(XTmp(:,jj), convkernel(kk,:))];  
    end
end
oldXTmp = XTmp;
XTmp = [XTmp, XConv(31:length(XTmp)+30,:)];


%% add airpuff Triggers

PuffT = zeros( length(Vidtime.(sess).wcam), 1);
for jj=1:p.nTrials
    t1 = find(Vidtime.(sess).wcam(:,2)>Puff.(sess)(jj,1),1);
    t2 = find(Vidtime.(sess).wcam(:,2)>Puff.(sess)(jj,2),1);
    PuffT(t1:t2-1)=1;
end
for jj=1:9
    PuffConv(:,jj) = conv(PuffT, convkernel2(jj,:)); 
end

XTmp=[PuffConv(1:length(XTmp), :), XTmp];

%% resample
XTimes = Vidtime.(sess).wcam(:,2);
[xtm, ind] = unique( XTimes);

dFFTimes = Norm.(sess).t;
t_int = [dFFTimes(1) dFFTimes(end,1)];

newX = [];
for jj=1:size(XTmp, 2)
   [tmp, newTimes] = resample_ds( XTmp(ind,jj), xtm, ds_rate, t_int ); 
   newX = [newX, tmp'];
end

[loco, ~] = resample_ds( -smooth(Speed.(sess)(:,2), 50), Speed.(sess)(:,1), ds_rate, t_int);
newX = [newX, 0.2*zscore(loco)'];

%% Zscore all except puff
newX = [ones(length(newX),1), newX];
newX(newX<1e-3)=0;


newDff = [];
newDff_noPC1 = [];
pccomp = [];

for jj=1:p.nROIs
    tmpNorm(:,jj) = inpaint_nans(Norm.(sess).n(:,jj)); 
end
noPC1 = subspace_svd(tmpNorm', -1 )';
[~,allpc] = pca(tmpNorm);

for jj=1:p.nROIs
   [tmp, ~] = resample_ds( tmpNorm(:,jj), Norm.(sess).t(:,jj), ds_rate, t_int );
    newDff = [newDff, tmp'];
    
   [tmp, ~] = resample_ds( noPC1(:,jj), Norm.(sess).t(:,jj), ds_rate, t_int );
    newDff_noPC1 = [newDff_noPC1, tmp'];

    [tmp,~] = resample_ds( allpc(:,jj), nanmean(Norm.(sess).t,2), ds_rate, t_int);
    pccomp = [pccomp, tmp'];
end

newDff( Seg.(sess).error, :)= NaN;
newDff_noPC1( Seg.(sess).error, :)= NaN;

newDff( error_frames, :)= NaN;
newDff_noPC1( error_frames, :)= NaN;
pccomp( error_frames,:) =NaN;

newX(error_frames, :)=NaN;

%% regression
clear b b2 b3
for jj=1:p.nROIs
   [b.All(:,jj), ~, Res(:,jj)] = regress( newDff(:, jj), newX ); 
   [b.noPC(:,jj), ~, Res2(:,jj)] = regress( newDff_noPC1(:, jj), newX ); 
   [b.PC(:,jj), ~, Res3(:,jj)] = regress( pccomp(:, jj), newX ); 
   
end
Y.All = newDff;
Y.noPC = newDff_noPC1;
Y.PC = pccomp;


%% compute residuals per variable
type = 'PC';

allcat = {[2:10], [11:13,17:43], [14:15,44:61],[16,62:70], 71};
allcatnames ={'Puff', 'WAngle', 'WPad', 'WMI', 'Speed'};
nCat = length(allcat);
nFactors = size(newX,2);

unqRS = nan(nCat, p.nROIs);
unqadjRS = unqRS;
allRS = nan(1, p.nROIs);
alladjRS = nan(1,p.nROIs);

for roi=1:p.nROIs
   ypred = newX*b.(type)(:,roi);
   sse = nansum((Y.(type)(:, roi)-ypred).^2);
   sst = nansum((Y.(type)(:, roi)-nanmean(Y.(type)(:, roi))).^2);
   nTm = sum(~isnan(Y.(type)(:,roi)));
   allRS(roi) = 1-sse/sst;
   alladjRS(roi) = 1-(nTm-1)/(nTm-nFactors-1)*sse/sst;
   for jj = 1:nCat
       use_ids = true(1,nFactors); use_ids(allcat{jj}) = false;
       ypred = newX(:,use_ids)*b.(type)(use_ids,roi);
       sse = nansum((Y.(type)(:, roi)-ypred).^2);
       unqRS(jj,roi) = (allRS(roi) - (1-sse/sst));
       unqadjRS(jj,roi) = (alladjRS(roi) - (1- (nTm-1)/(nTm-nFactors-2)*sse/sst));
   end
    
end
unqRS = unqRS./sum(unqRS,1);
unqadjRS = unqadjRS./sum(unqadjRS,1);
