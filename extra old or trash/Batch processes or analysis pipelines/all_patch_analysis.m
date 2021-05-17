% %% get raw data
% 
% [params, raw, times, ~] = load_POI();
% 
% %-------------------------------------------------------------------------%
%% Do Motion correction
% Fixed Ref
%-----------------------
% montage = false;
% nX = size(raw, 2);
% nY = size(raw,1);
% Tm = size(raw,3);
% nTrials = size(raw,4);
% nROI = size(raw, 5);
% 
% DX =4; DY = 4;
% Fx = DX+1:nX-DX;
% Fy = DY+1:nY-DY;
% 
% nx = size(Fx,2);
% ny = size(Fy,2);
% 
% [ all_disp_traj, refs, MC_correct ] = complete_mc( raw, [DX, DY], 0.2, 1 );
% disp_x = ceil( all_disp_traj/(2*DY+1) ) - DX-1;
% disp_y = mod( all_disp_traj-1, 2*DY+1 ) - DY;

% Slow Ref
%-----------------------
% montage = false;
% nX = size(raw, 2);
% nY = size(raw,1);
% Tm = size(raw,3);
% nTrials = size(raw,4);
% nROI = size(raw, 5);
% 
% DX =4; DY = 4;
% Fx = DX+1:nX-DX;
% Fy = DY+1:nY-DY;
% 
% nx = size(Fx,2);
% ny = size(Fy,2);
% 
% [ all_disp_traj, refs, MC_correct ] = complete_mc_slow_ref( raw2, [DX, DY], 0.2, 1 );
% disp_x = ceil( all_disp_traj/(2*DY+1) ) - DX-1;
% disp_y = mod( all_disp_traj-1, 2*DY+1 ) - DY;

% With Montage
% -----------------------
montage = true;
nX = size(raw, 2);
nY = size(raw,1);
Tm = size(raw,3);
nTrials = size(raw,4);
nROI = size(raw, 5);

DX =4; DY = 4;
Fx = DX+1:nX-DX;
Fy = DY+1:nY-DY;

nx = size(Fx,2);
ny = size(Fy,2);

[ all_disp_traj, refs, MC_correct ] = complete_mc( raw, [DX, DY], 0.2, 1 , 'gamma', 1, 'montage', true, 'change_ref', true);
disp_x = ceil( all_disp_traj/(2*DY+1) ) - DX-1;
disp_y = mod( all_disp_traj-1, 2*DY+1 ) - DY;

%-------------------------------------------------------------------------%
%% get object masks for original and MC-stack

disp('Creating object masks' )
im={};
for roi=1:nROI
    im{roi} = mean( reshape(MC_correct(:,:,:,:,roi), [ny,nx, Tm*nTrials]), 3 );
end

masks={};
for jj=1:nROI
[~,threshold]=edge(im{jj},'sobel');
BW=edge(im{jj}, 'sobel',0.4*threshold);
BWdil = imdilate(BW, strel('square',3));
masks{jj} = imfill(BWdil, 'holes');
end


im2={};
for roi=1:nROI
    im2{roi} = mean( reshape(raw(Fy,Fx,:,:,roi), [ny,nx, Tm*nTrials]), 3 );
end

masks2={};
for jj=1:nROI
[~,threshold]=edge(im2{jj},'sobel');
BW=edge(im2{jj}, 'sobel',0.4*threshold);
BWdil = imdilate(BW, strel('square',3));
masks2{jj} = imfill(BWdil, 'holes');
end


%-------------------------------------------------------------------------%
%% apply masks

disp('Applying masks to patches...' )
MC_and_mask = nan(Tm, nTrials, nROI);
time_mc_and_mask = MC_and_mask;

noMC_and_mask = MC_and_mask;
time_nomc_and_mask = MC_and_mask;

% For MC
parfor roi=1:nROI
   temp=reshape(MC_correct(:,:,:,:,roi), [ny*nx, Tm, nTrials]);
   
   for trial = 1:nTrials
   for t=1:Tm
       mc = (roi + (t-1)*nROI)*(1-montage) + montage*t;
       MC_and_mask(t,trial,roi) = mean( temp(masks{jj}(:), t, trial), 1 );
       timetemp = squeeze(times(Fy+disp_y(trial, mc), Fx+disp_x(trial, mc), t, trial, roi));
       time_mc_and_mask(t, trial, roi) = mean( timetemp(masks{jj}) );
   end
   end
end

% For no MC
parfor roi=1:nROI
   rawtemp = reshape(raw(Fy, Fx,:,:,roi), [ny*nx, Tm, nTrials]);
   timetemp = reshape(times(Fy, Fx,:,:,roi), [ny*nx, Tm, nTrials]);
   for trial = 1:nTrials
   for t=1:Tm
       noMC_and_mask(t,trial,roi) = mean( rawtemp(masks{jj}(:), t, trial) );
       time_nomc_and_mask(t, trial, roi) = mean( timetemp(masks{jj}(:), t, trial) );
   end
   end
end


%-------------------------------------------------------------------------%
%% Normalize and smoothen
params.baseline_per_trial=false;
params.norm_method='percentile';
params.base_percentile=20;
params.smooth_type='Moving average';
params.smooth_scale=60;
disp('Normalisation and smoothening...' )
[norm, F0] =normalize_and_smooth( MC_and_mask, time_mc_and_mask, params);


%-------------------------------------------------------------------------%
%% Get speed data and locomotion periods
disp('Speed data and locomotion periods...' )
fp = [params.exp_path, '/Speed_Data/'];
speed = get_speed_data( fp);
spd = [];
for jj=1:nTrials
    spd = [spd; speed{jj}];
end

lp= get_loco_period( speed, true, 50, 3000, [0.5 0.2], 30)

%-------------------------------------------------------------------------%
%% Plot traces and locomotion
figure;
cmap=colormap(lines);
for jj=1:size(lp,1)
hold on
a=area(lp(jj,1:2),repmat(8+nROI,[1,2]));
a.FaceAlpha=0.25;
a.EdgeAlpha=0;
end
plot_traces( norm, time_mc_and_mask, [1 3 2])
hold on
for jj=1:nTrials
plot( speed{jj}(:,1), 0.5*smooth( abs(speed{jj}(:,2)), 100 )+4+nROI)
end


%-------------------------------------------------------------------------%
%% Lag correlations between locomotion and dFF
disp('Calculating correlations between speed and dFF during locomotion periods...' )
[lags,corrs,~]=all_period_lagcorr(lp, abs(spd(:,2)), spd(:,1), reshape( norm, [Tm*nTrials,nROI]), reshape(time_mc_and_mask, [Tm*nTrials,nROI]),50,1000,500);

%-------------------------------------------------------------------------%
%% Plot masks and Mean MC-image

ncols=4;
nrows=5;

%------------------
% With MC
for jj=1:nROI
if mod(jj,ncols*nrows)==1
        figure()
        suptitle(sprintf('Mean Image(with MC) - Patches %d to %d', jj, min(jj-1+nrows*ncols, nROI)))
end
subplot(nrows,ncols,mod(jj-1,nrows*ncols)+1)
imagesc(im{jj})
end   


for jj=1:nROI
if mod(jj,ncols*nrows)==1
        figure()
        suptitle(sprintf('Masks(with MC) - Patches %d to %d', jj, min(jj-1+nrows*ncols, nROI)))
end
subplot(nrows,ncols,mod(jj-1,nrows*ncols)+1)
imagesc(masks{jj})
end   

%------------------
% Without MC
for jj=1:nROI
if mod(jj,ncols*nrows)==1
        figure()
        suptitle(sprintf('Mean Image(without MC) - Patches %d to %d', jj,min(jj-1+nrows*ncols, nROI)))
end
subplot(nrows,ncols,mod(jj-1,nrows*ncols)+1)
imagesc(im2{jj})
end   


for jj=1:nROI
if mod(jj,ncols*nrows)==1
        figure()
        suptitle(sprintf('Masks(without MC) - Patches %d to %d', jj, min(jj-1+nrows*ncols, nROI)))
end
subplot(nrows,ncols,mod(jj-1,nrows*ncols)+1)
imagesc(masks2{jj})
end   
   
   