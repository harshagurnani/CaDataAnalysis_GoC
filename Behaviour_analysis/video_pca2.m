% close all
% clear variables

clear vid vid_motion video video_motion

%%
file='forepawR_ecam';
width = 358; height = 218;

folder   = 'D:\myStuff\data_copy_from_server\HG 03 Caroline Herschel\exp 03\videos\180423_12_07_09\split';
save_dir = 'D:\myStuff\data_copy_from_server\HG 03 Caroline Herschel\exp 03\videos\180423_12_07_09\';


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

if ~exist( 'find_eig', 'var' ), find_eig = true; end

vid         = double(reshape(vid, [width*height, nframes]));
vid_motion  = diff(vid,1,2);

if find_eig
    [coeff,     score,      latent]     = pca(vid,          'NumComponents',pca_numcomp);
    [coeff_m,   score_m,    latent_m]   = pca(vid_motion,   'NumComponents',pca_numcomp);

    video.comp          = score;    video.lowd_data         = coeff';    video.var           = latent;
    video_motion.comp   = score_m;  video_motion.lowd_data  = coeff_m';  video_motion.var    = latent_m;

else
    vid         = vid - mean(vid,2);    %Centering
    vid_motion  = vid_motion - mean(vid_motion,2);

    % find normalisation factors
    pca_numcomp = size(eigenvec,2);
    if ~exist( 'eig_var', 'var')
        eig_var     = arrayfun(@(col) norm(eigenvec(:,col)),    1:pca_numcomp);
        eig_var_m   = arrayfun(@(col) norm(eigenvec_m(:,col)),  1:pca_numcomp);
    end
    video.comp      = eigenvec;    %should exist in workspace
    video.lowd_data = eigenvec'*vid./repmat((eig_var.^2)',1,nframes);
    video.var       = latent;
    
    video_motion.comp      = eigenvec_m;
    video_motion.lowd_data = eigenvec_m'*vid_motion./repmat((eig_var_m.^2)',1,nframes-1);
    video_motion.var       = latent_m;
end
%% Saving

save( [save_dir, file, '_pca.mat'], 'video', 'video_motion', '-v7.3')

%% view results

fp='forepawR ecam';
width = 358; height = 218;
pca_numcomp = 100;

score = reshape(video.comp, [height, width, pca_numcomp]);
score_m = reshape(video_motion.comp, [height, width, pca_numcomp]);

figure; for jj=1:25, subplot(5,5,jj); imagesc(score(:,:,jj), [-100, 100]); end
suptitle( ['Components of video from ', fp])

figure; for jj=1:25, subplot(5,5,jj); imagesc(score_m(:,:,jj), [-100, 100]); end
suptitle( ['Components of motion from ', fp])

coeff = video.lowd_data; 
figure; cm=colormap(jet(10)); for jj=1:10, plot( coeff(jj,:) + 0.2*jj,'color', cm(jj,:) ); hold on; end
title(['Evolution of top 10 components of video from ', fp])

coeff_m = video_motion.lowd_data;
figure; cm=colormap(jet(10)); for jj=1:10, plot( coeff_m(jj,:) + 0.2*jj,'color', cm(jj,:) ); hold on; end
title(['Evolution of top 10 components of motion from ', fp])

