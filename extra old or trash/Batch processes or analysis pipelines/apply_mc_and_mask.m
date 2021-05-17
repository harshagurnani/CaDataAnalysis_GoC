im = nan(52,29,Tm, nTrials, nROI); 
roin = [1 2 3 4 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21];
[xyshift, rho, bestimage ] = MC_with_ECC( raw(:,:,:,:,roin), ref, 'levels', 2);
xcrop = 3:31;
ycrop = 4:55;
nx = size(xcrop,2);
ny = size(ycrop,2);

nROI = size(roin,2);
for roi=1:nROI
for trial=1:nTrials
for t=1:Tm
[wim, ~] = iat_inverse_warping(raw(:,:,t,trial,roin(roi)), xyshift(:,t, trial), par.transform, xcrop, ycrop, 'linear');
im(:,:,t,trial,roi) = wim;
end
end
end
MC_correct = im;


im={};
for roi=1:nROI
im{roi} = mean( reshape(MC_correct(:,:,:,:,roi), [52,29, Tm*nTrials]), 3 );
end
masks={};
for jj=1:nROI
[~,threshold]=edge(im{jj},'sobel');
BW=edge(im{jj}, 'sobel',0.5*threshold);
BWdil = imdilate(BW, strel('square',3));
masks{jj} = imfill(BWdil, 'holes');
end

for roi=1:nROI
temp=reshape(MC_correct(:,:,:,:,roi), [52*29, Tm, nTrials]);
for trial=1:nTrials
for t=1:Tm
MC_and_mask(t,trial,roi) = mean( temp(masks{jj}(:), t, trial), 1 );
end
end
end
time_mn = mean(reshape(times(:,:,:,:,roin), [58*33,Tm,nTrials, nROI]),1);
time_mn=reshape(time_mn, [Tm, nTrials, nROI]);

params.baseline_per_trial=false;
params.norm_method='percentile';
params.base_percentile=5;
params.smooth_type='Moving average';
params.smooth_scale=100;
[norm, F0] =normalize_and_smooth( MC_and_mask, time_mn, params);

plot_traces(norm, time_mn, [1 3 2])