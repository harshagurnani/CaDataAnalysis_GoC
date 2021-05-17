%% Baseline
% Load data
[p,r,t,~] = load_POI('source',exp.baseline,     'ROI_type','patch',     'avg_patch',false,  'normalize', false );
[nL, nP, nTm, nTr, nROI ] = size(r);
r = reshape(r, [nL, nP, nTm*nTr, nROI]);
t=reshape( t, [nL, nP, nTm*nTr, nROI] );
clear MC

% select ref and do Motion Correction
figure; nr = ceil(sqrt(nROI));
for roi=1:nROI
subplot(nr,nr,roi);     imagesc( nanmean(r(:,:,:,roi), 3 ), [50 500]);     title(num2str(roi));
end

MC.ref_ROI  = input('Enter ROI numbers (as vector) to use for motion correction. Return empty if no need to correct:   ');
% MC.ref_ROI=[1 2 3 5 7 10 22];
MC.max_disp = 8;

if ~isempty(MC.ref_ROI)
    if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0
        MC.ref_ROI = 1:nROI;
        [MC.ref, MC.baseline.offsets, r] = do_mc_fastcorr( r, 'do_crop', true, 'max_disp', MC.max_disp   );
        [~,~,t]                          = do_mc_fastcorr( t, 'do_crop', true, 'use_offsets', MC.baseline.offsets );
    else
        size_needed = length(MC.ref_ROI)*nL*nP*nTm*nTr*16 / (1024*1024*1024);  % gigabytes
        if size_needed < 0.5
            [MC.ref, MC.baseline.offsets   ] = do_mc_fastcorr( r(:,:,:,MC.ref_ROI) );
            for roi=1:nROI
                [ ~,~, r(:,:,:,roi)]               = do_mc_fastcorr( r(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets, 'max_disp', MC.max_disp );
                [ ~,~, t(:,:,:,roi)]               = do_mc_fastcorr( t(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets, 'max_disp', MC.max_disp );
            end
        else
            % do per trial - too much data to load simultaneously
            r = reshape(r, [nL, nP, nTm, nTr, nROI]);
            t = reshape(t, [nL, nP, nTm, nTr, nROI]);
            MC.baseline.offsets = cell(nTr,1);
            [MC.ref, MC.baseline.offsets{1}   ] = do_mc_fastcorr( r(:,:,:,1,MC.ref_ROI) );
            for roi=1:nROI
                [ ~,~, r(:,:,:,1,roi)]               = do_mc_fastcorr( r(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets{1}, 'max_disp', MC.max_disp );
                [ ~,~, t(:,:,:,1,roi)]               = do_mc_fastcorr( t(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets{1}, 'max_disp', MC.max_disp );
            end

            for trial = 2:nTr
                [~, MC.baseline.offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                for roi=1:nROI
                    [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets{trial}, 'max_disp', MC.max_disp );
                    [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.baseline.offsets{trial}, 'max_disp', MC.max_disp );
                end
                
            end
            r = reshape(r, [nL, nP, nTm*nTr, nROI]);
            t=reshape( t, [nL, nP, nTm*nTr, nROI] );
            MC.baseline.offsets = cell2mat(MC.baseline.offsets);
        end
    end
end
% show result
%=
MC.baseline.error = ( abs( MC.baseline.offsets(:,1) ) > MC.max_disp | abs( MC.baseline.offsets(:,2) ) > MC.max_disp );

Basic.baseline = struct( 'r',r,'t',t,'params',p );
save( [exp.analysed, '\raw_', suffix.baseline, '.mat'], 'Basic', 'MC', 'exp', '-v7.3' );
clear Basic

disp('.....................Saved Baseline Data after motion correction')
MC_show = strcmpi(input('Do you want to test MC quality? y/n:    ','s'), 'y');
if MC_show
    for tm=sort(unidrnd(nTm*nTr, 1,50))
    for roi=1:nROI
        subplot(nr,nr,roi);  imagesc(r(:,:,tm,roi),[50 2000]);   
    end
    suptitle(['Frame ',num2str(tm)]);
    pause(0.002)
    end
end

%% ----------------------------------------------------------------------
if do_AP
% load AP experiment
[p2,r2,t2,~] = load_POI('source',exp.AP,     'ROI_type','patch',     'avg_patch',false,  'normalize', false );
[nL, nP, nTm, nTr, nROI ] = size(r2);
r2 = reshape(r2, [nL, nP, nTm*nTr, nROI]);
t2 = reshape(t2, [nL, nP, nTm*nTr, nROI]);

% Align to same reference
if ~isempty(MC.ref_ROI)
    if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0
        MC.ref_ROI = 1:nROI;
        [~, MC.AP.offsets, r2] = do_mc_fastcorr( r2, 'do_crop', true, 'ref', MC.ref, 'max_disp', MC.max_disp );
        [~,~,t2]                    = do_mc_fastcorr( t2, 'do_crop', true, 'use_offsets', MC.AP.offsets, 'max_disp', MC.max_disp );
    else
        if size_needed < 0.5
            [~, MC.AP.offsets   ] = do_mc_fastcorr( r2(:,:,:,MC.ref_ROI), 'ref', MC.ref );
            for roi=1:nROI
                [ ~,~, r2(:,:,:,roi)]                 = do_mc_fastcorr( r2(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.AP.offsets, 'max_disp', MC.max_disp );
                [ ~,~, t2(:,:,:,roi)]                 = do_mc_fastcorr( t2(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.AP.offsets, 'max_disp', MC.max_disp );
            end
        else
            r2 = reshape(r2, [nL, nP, nTm, nTr, nROI]);
            t2 = reshape(t2, [nL, nP, nTm, nTr, nROI]);
            MC.AP.offsets = cell(nTr,1);
            
            for trial = 1:nTr
                [~, MC.AP.offsets{trial}   ] = do_mc_fastcorr( r2(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                for roi=1:nROI
                    [ ~,~, r2(:,:,:,trial,roi)]               = do_mc_fastcorr( r2(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.AP.offsets{trial}, 'max_disp', MC.max_disp );
                    [ ~,~, t2(:,:,:,trial,roi)]               = do_mc_fastcorr( t2(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.AP.offsets{trial}, 'max_disp', MC.max_disp );
                end
                
            end
            r2 = reshape(r2, [nL, nP, nTm*nTr, nROI]);
            t2 = reshape(t2, [nL, nP, nTm*nTr, nROI] );
            MC.AP.offsets = cell2mat(MC.AP.offsets);
        end   
    end
end
MC.AP.error = ( abs( MC.AP.offsets(:,1) ) > MC.max_disp | abs( MC.AP.offsets(:,2) ) > MC.max_disp );

Basic.baseline = struct( 'r',r2,'t',t2,'params',p2 );
save( [exp.analysed, '\raw_', suffix.AP, '.mat'], 'Basic', 'MC','exp' , '-v7.3' );
clear Basic MC_show
disp('.....................Saved Air Puff Imaging Data after motion correction')
end

%% Masks and Segmentation
if do_AP
    joint_data = cat(3, r, r2);
else
    joint_data = r;
end
im      = arrayfun(@(roi) nanmean(joint_data(:,:,:,roi),3),1:nROI,'UniformOutput',false);           % mean image across both experiments
keeppix = arrayfun(@(roi) ~isnan(mean(joint_data(:,:,:,roi),3)),1:nROI,'UniformOutput',false);      % pixels present in both experiments

Seg.method = 'sobel';
clear joint_data
for roi=1:nROI
[~,threshold]=edge(im{roi}, Seg.method);    Seg.threshold(roi) = 0.5*threshold;
BW      = edge(im{roi}, Seg.method, Seg.threshold(roi) );
BWdil   = imdilate(BW, strel('square',3));
Seg.fullmasks{roi} = imfill(BWdil, 'holes').*keeppix{roi};

% keep maximal connected components
Seg.masks{roi} = false(size(Seg.fullmasks{roi}));
Seg.masks{roi}( get_bigg_conncomp(Seg.fullmasks{roi}) ) = true;
end
Seg.mean_image = im;
clear keeppix BW BWdil

[ Seg.baseline.r,  Seg.baseline.t ] = mask_nonbinary_and_avg( r, t, Seg.masks);
% remove uncorrected frames
errorid = find(MC.baseline.error);
Seg.baseline.error = MC.baseline.error;
if max(diff(errorid))<8, Seg.baseline.error(errorid(1):errorid(end)) = true; end
Seg.baseline.r( Seg.baseline.error,: ) = NaN;

if do_AP
    [ Seg.AP.r,        Seg.AP.t       ] = mask_nonbinary_and_avg( r2, t2, Seg.masks);
    errorid = find(MC.AP.error);
    Seg.AP.error = MC.AP.error;
    if max(diff(errorid))<8, Seg.AP.error(errorid(1):errorid(end))=true; end
    Seg.AP.r( Seg.AP.error,: ) = NaN;
    clear r2 t2
end
clear r t errorid threshold nL nP nTm nTr 

save( [exp.analysed, '\masked_raw.mat'], 'Seg', 'exp', '-v7.3' );
disp('.....................Saved Masked Data')

%%
% Seg.method    =intensity threshold';
% Seg.threshold = 250;
% Seg.masks = arrayfun( @(roi) mean(joint_data(:,:,:,roi),3) > Seg.threshold, 1:nROI, 'UniformOutput', false );

%% Normalisation
% in paint NaNs only for normalization process
for roi=1:nROI
    Seg.baseline.r(:,roi) = inpaint_nans( Seg.baseline.r(:,roi),4);  
end
if do_AP
for roi=1:nROI
    Seg.AP.r(:,roi)       = inpaint_nans( Seg.AP.r(:,roi),4);  
end
end

% use defaults for now
p.smooth_causal = true; if do_AP, p2.smooth_causal = true; end

Norm.baseline.n   = normalize_and_smooth( Seg.baseline.r, Seg.baseline.t, p );
Norm.baseline.t   = Seg.baseline.t; 
Norm.baseline.n(Seg.baseline.error,:) = NaN;
MC_error.baseline = MC.baseline.error;

if do_AP
Norm.AP.n         = normalize_and_smooth( Seg.AP.r, Seg.AP.t, p2 );
Norm.AP.t = Seg.AP.t;
Norm.AP.n(Seg.AP.error,:) = NaN;
MC_error.AP = MC.AP.error;
else
    p2=[];
end


masks = Seg.masks;
save([exp.analysed, '\masked_norm_test.mat'], 'Norm','p','p2','masks','im','MC_error','-v7.3');

clear MC_error masks im

disp('.....................Saved Normalised Traces')

%% HELPER
%------------------------------------------------------------------------%
function pix = get_bigg_conncomp( bin_img )
components     = bwconncomp(bin_img);
numPixels      = cellfun(@numel,components.PixelIdxList);
[~,idx]        = max(numPixels);
pix            = components.PixelIdxList{idx(1)};
end

