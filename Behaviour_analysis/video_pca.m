% close all
% clear variables
folder = 'D:\myStuff\data_copy_from_server\HG 03 Caroline Herschel\exp 03\videos\180423_13_22_39\split';
file='forepawR_ecam';
width = 358; height = 218;
save_dir = 'D:\myStuff\data_copy_from_server\HG 03 Caroline Herschel\exp 03\videos\180423_13_22_39\';


%% Load all files
allfiles    = dir([folder,'\*',file,'*']);
nfiles      = numel(allfiles);

vid = cell(1,1,nfiles);
for fn=1:nfiles
    tmp = imfinfo([folder,'\',allfiles(fn).name] );
    nframes = numel(tmp); nW = tmp(1).Width; nH = tmp(1).Height;
    vid{fn} = uint8(nan(nH, nW, nframes)); 
    for frame = 1:nframes
    vid{fn}(:,:,frame) = imread( [folder,'\',allfiles(fn).name], frame );
        
    end
end
vid = cell2mat(vid);
if nW~=width || nH~=height
    vid = vid(1:height, 1:width, :);
end
clear allfiles tmp nW nH nfiles frame fn 

%% PCA
tmp = size(vid); nframes = tmp(end); clear tmp
pca_numcomp = 100;

vid         = double(reshape(vid, [width*height, nframes]));
vid_motion  = diff(vid,1,2);
[coeff,     score,      latent]     = pca(vid,          'NumComponents',pca_numcomp);
[coeff_m,   score_m,    latent_m]   = pca(vid_motion,   'NumComponents',pca_numcomp);

%% Saving
video.comp          = score;    video.lowd_data         = coeff;    video.var           = latent;
video_motion.comp   = score_m;  video_motion.lowd_data  = coeff_m;  video_motion.var    = latent_m;

save( [save_dir, file, '_pca.mat'], 'video', 'video_motion', '-v7.3')

%% view results
score = reshape(video.comp, [height, width, pca_numcomp]);
score_m = reshape(video_motion.comp, [height, width, pca_numcomp]);
figure; for jj=1:25, subplot(5,5,jj); imagesc(score(:,:,jj), [-100, 100]); end
suptitle( 'Components of video')

figure; for jj=1:25, subplot(5,5,jj); imagesc(score_m(:,:,jj), [-100, 100]); end
suptitle( 'Components of motion')

coeff = video.lowd_data; 
figure; cm=colormap(jet(10)); for jj=1:10, plot( coeff(:,jj) + 0.2*jj,'color', cm(jj,:) ); hold on; end
title('Evolution of top 10 components of video')

coeff_m = video_motion.lowd_data;
figure; cm=colormap(jet(10)); for jj=1:10, plot( coeff_m(:,jj) + 0.2*jj,'color', cm(jj,:) ); hold on; end
title('Evolution of top 10 components of motion')

