do_ref  = false;
ask_plot = false;%true;%;

%% setup ref for movement correction - if already specified
if ~isempty( exp.Folder.UseRef )
    load( exp.Folder.UseRef, 'MC' );
    ref = MC.ref;   ref_roi = MC.ref_ROI; do_ref = false;
    clear MC
    MC = struct( 'ref', ref, 'ref_ROI', ref_ROI' );
end

%% For all sessions - Load and do Motion Correction
for sn=1:nSess

sess            = session_types{use_sessions(sn)};
suffix          = sess;

% Load data
    [p,r,t,~] = load_POI('source',exp.(sess),     'ROI_type','patch',     'avg_patch',false,  'normalize', false );
    [nL, nP, nTm, nTr, nROI ] = size(r);
    r = reshape(r, [nL, nP, nTm*nTr, nROI]);
    t=reshape( t, [nL, nP, nTm*nTr, nROI] );

MC.max_disp = 8;

if do_ref
% select ref and do motion correction
    f=figure; nr = ceil(sqrt(nROI));
    for roi=1:nROI
    subplot(nr,nr,roi);     imagesc( nanmean(r(:,:,:,roi), 3 ));     title(num2str(roi)); %, [50 500]
    end
    MC.ref_ROI  = input('Enter ROI numbers (as vector) to use for motion correction. Return empty if no need to correct:   ');
    close(f)
    %%%%%%%%%%%%%%%%%%%%5
    
    if ~isempty(MC.ref_ROI)
        if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0
            MC.ref_ROI = 1:nROI;
            % one shot reference and MC
            [MC.ref, MC.(sess).offsets, r] = do_mc_fastcorr( r, 'do_crop', true, 'max_disp', MC.max_disp, 'speed', Speed.(sess), 'times', t(1,1,:,1)   );
            [~,~,t]                        = do_mc_fastcorr( t, 'do_crop', true, 'use_offsets', MC.(sess).offsets );
        else
            size_needed = length(MC.ref_ROI)*nL*nP*nTm*nTr*16 / (1024*1024*1024);  % gigabytes
            if size_needed < 0.2
                %creare ref and get offsets
                [MC.ref, MC.(sess).offsets   ] = do_mc_fastcorr( r(:,:,:,MC.ref_ROI), 'max_disp', MC.max_disp, 'speed', Speed.(sess), 'times', t(1,1,:,1)  ); 
                %apply offsets
                for roi=1:nROI
                    [ ~,~, r(:,:,:,roi)]               = do_mc_fastcorr( r(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                    [ ~,~, t(:,:,:,roi)]               = do_mc_fastcorr( t(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                end
            else
                % do per trial - too much data to load simultaneously
                r = reshape(r, [nL, nP, nTm, nTr, nROI]);
                t = reshape(t, [nL, nP, nTm, nTr, nROI]);
                MC.(sess).offsets = cell(nTr,1);
                % create ref using noloco /trial 1
                trial = 1;
                %% ind trial with no loco
                speed_dt = (Speed.(sess)(2,1)-Speed.(sess)(1,1));  bin5s = min( floor(1500/speed_dt), length(Speed.(sess)) ); 
                sumspd = arrayfun( @(jj) sum(abs(Speed.(sess)(jj+1:jj+bin5s,2) )), 0:length(Speed.(sess))-bin5s );    
                noloco = find(sumspd == 0, 1); noloco_time = noloco*speed_dt; 
                TrialTimes = squeeze(t(1,1,1,:,1));
                trial = find(noloco_time<TrialTimes,1)-1; if isempty(trial), trial = nTr; end
                t1 = find(Speed.(sess)(:,1)>TrialTimes(trial),1); t2 =  find(Speed.(sess)(:,1)>TrialTimes(trial+1),1); if isempty(t2), t2 = length(Speed.(sess)); end
%                 [MC.ref, MC.(sess).offsets{1}   ] = do_mc_fastcorr( r(:,:,:,1,MC.ref_ROI), 'max_disp', MC.max_disp, 'speed', Speed.(sess), 'times', t(1,1,:, trial,1)  );
                [MC.ref, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'max_disp', MC.max_disp, 'speed', Speed.(sess)(t1:t2-1,:), 'times', t(1,1,:, trial,1)  );
                % apply offsets on trial 1
%                 for roi=1:nROI
%                     [ ~,~, r(:,:,:,1,roi)]               = do_mc_fastcorr( r(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
%                     [ ~,~, t(:,:,:,1,roi)]               = do_mc_fastcorr( t(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
%                 end

                for trial = 1:nTr %instead of 2:nTr
                    % get offsets
                    [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                    % apply offsets
                    for roi=1:nROI
                        [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp, 'newreg_default', 'nans' );
                        [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                    end

                end
                r = reshape(r, [nL, nP, nTm*nTr, nROI]);
                t=reshape( t, [nL, nP, nTm*nTr, nROI] );
                MC.(sess).offsets = cell2mat(MC.(sess).offsets);
            end
        end
    else
        MC.(sess).offsets = zeros( nTm*nTr,2);
        MC.ref = [];
    end
%     ask_plot = true;
    do_ref = false;
    
else
% Align to same reference
    if ~isempty(MC.ref_ROI)
        if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0
            MC.ref_ROI = 1:nROI;
            %get and apply offsets
            [~, MC.(sess).offsets, r]  = do_mc_fastcorr( r, 'do_crop', true, 'ref', MC.ref, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
            [~,~,t]                    = do_mc_fastcorr( t, 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
        else
            if size_needed < 0.2
                % get offsets
                [~, MC.(sess).offsets   ] = do_mc_fastcorr( r(:,:,:,MC.ref_ROI), 'ref', MC.ref );
                % apply offsets
                for roi=1:nROI
                    [ ~,~, r(:,:,:,roi)]                 = do_mc_fastcorr( r(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                    [ ~,~, t(:,:,:,roi)]                 = do_mc_fastcorr( t(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                end
            else
                r = reshape(r, [nL, nP, nTm, nTr, nROI]);
                t = reshape(t, [nL, nP, nTm, nTr, nROI]);
                
                MC.(sess).offsets = cell(nTr,1);

                for trial = 1:nTr
                    %get offsets
                    [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                    % apply offsets
                    for roi=1:nROI
                        [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp, 'newreg_default', 'nans'  );
                        [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                    end

                end
                r = reshape(r, [nL, nP, nTm*nTr, nROI]);
                t = reshape(t, [nL, nP, nTm*nTr, nROI] );
                MC.(sess).offsets = cell2mat(MC.(sess).offsets);
            end   
        end
        
    else
        MC.(sess).offsets = zeros( nTm*nTr,2);
    end
    
end

    
% check for error frames
MC.(sess).error = ( abs( MC.(sess).offsets(:,1) ) > MC.max_disp | abs( MC.(sess).offsets(:,2) ) > MC.max_disp );

% save data
Basic.(sess) = struct( 'r',r,'t',t,'params',p );
save( [exp.analysed, '/raw_', suffix, '.mat'], 'Basic', 'MC', 'exp', 'p', '-v7.3' );

% save mean image and pixels present
if ~exist('im','var')
    im = cell(nROI,1); for jj=1:nROI, im{jj} = nan(nL,nP,nSess); end
    keeppix = im;
end
for roi=1:nROI
    im{roi}(1:nL, 1:nP, sn)      = nanmean(r(:,:,:,roi),3);           % mean image
    keeppix{roi}(1:nL, 1:nP, sn) = ~isnan(mean(r(:,:,:,roi),3));      % pixels present in all frames
end

clear Basic
disp(['.....................Saved Imaging Data for ', sess, ' session after motion correction'])

% show result
if ask_plot
    MC_show = strcmpi(input('Do you want to test MC quality? y/n:    ','s'), 'y');
    if MC_show
        % show random N frames
        nframes = 50;
        for tm=sort(unidrnd(nTm*nTr, 1,nframes))
        for roi=1:nROI
            subplot(nr,nr,roi);  imagesc(r(:,:,tm,roi),[50 2000]);   
        end
        suptitle(['Frame ',num2str(tm)]);
        pause(0.002)
        end
    end
    ask_plot = false;
end

clear r t

end

%% ----------------------------------------------------------------------

%% Masks /Segmentation
clear MC
clear suffix

for roi=1:nROI
im{roi}      = nanmean(im{roi}, 3);
keeppix{roi} = (sum(keeppix{roi},3) == nSess) ;      % pixels present in all experiments
end
Seg.mean_image = im;


if isempty( exp.Folder.UseMasks)
    % Create masks
    Seg.method = 'sobel';
    for roi=1:nROI
    [~,threshold]=edge(im{roi}, Seg.method);    Seg.threshold(roi) = 0.5*threshold;
    BW      = edge(im{roi}, Seg.method, Seg.threshold(roi) );
    BWdil   = imdilate(BW, strel('square',3));
    Seg.fullmasks{roi} = imfill(BWdil, 'holes');%.*keeppix{roi};    Do this multiplication afterwards

    % keep maximal connected components
    Seg.masks{roi} = false(size(Seg.fullmasks{roi}));
    Seg.masks{roi}( get_bigg_conncomp(Seg.fullmasks{roi}) ) = true;
    Seg.masks{roi} = Seg.masks{roi} .*keeppix{roi};
    
    % Make background mask    
    Seg.all_background_masks{roi} =  ~(imdilate(Seg.fullmasks{roi}, strel('disk',4,4))); % inverse of detected ROI components
    Seg.all_background_masks{roi} = Seg.all_background_masks{roi} .*keeppix{roi};
    
    % Dark background
    bg_lim = prctile(Seg.mean_image{roi}(Seg.fullmasks{roi}),10);       % lowest 10th of fullmask
    Seg.background_masks{roi} = Seg.all_background_masks{roi} .*(Seg.mean_image{roi}<bg_lim);
    
    end   
    clear BW BWdil threshold
    
    
else
    % Use masks
    load( exp.Folder.UseMasks, 'Seg' );
    fullmasks = Seg.fullmasks; masks = Seg.masks;
    background = Seg.background_masks;
    clear Seg
    Seg = struct( 'fullmasks', fullmasks, 'masks', masks, 'background_masks', background );
    clear masks fullmasks background
end
 

%% Apply masks
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    suffix          = sess;
    load( [exp.analysed, '/raw_', suffix, '.mat'], 'Basic', 'MC' );

    [ Seg.(sess).r,  Seg.(sess).t ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.masks);           % avg somatic fluorescence
    [ Seg.(sess).background, ~ ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.background_masks);   % avg dark baclground fluorescence
    
    % remove uncorrected frames
    errorid = find(MC.(sess).error);
    Seg.(sess).error = MC.(sess).error;
    if max(diff(errorid))<8, Seg.(sess).error(errorid(1):errorid(end)) = true; end
    Seg.(sess).r( Seg.(sess).error,: ) = NaN;
    
    clear Basic MC
end

save( [exp.analysed, '/masked_raw.mat'], 'Seg', 'exp', 'keeppix', '-v7.3' );
disp('.....................Saved Masked Data')
clear errorid keeppix

%% Normalisation
% in paint NaNs only for normalization process

for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    suffix          = sess;
    for roi=1:nROI
        Seg.(sess).r(:,roi) = inpaint_nans( Seg.(sess).r(:,roi),4);  
    end

    % use defaults for now (100 ms smoothing with moving average)
    load( [exp.analysed, '/raw_', suffix, '.mat'], 'p' );
    
    p.smooth_causal = false;
    if ~isfield( exp, 'params'), p.smooth_scale = 100;%ms %default
    else,                        p.smooth_scale  = exp.params.dff_smooth_scale; end
    nTm = p.nTimepoints; nTr = size(Seg.(sess).t,1)/nTm;
    
    % normalized soma mask
    [Norm.(sess).n, Norm.(sess).F0]   = normalize_and_smooth( reshape(Seg.(sess).r,[nTm,nTr,nROI]), reshape(Seg.(sess).t,[nTm,nTr,nROI]), p );
    Norm.(sess).n   = reshape(Norm.(sess).n, [nTm*nTr,nROI]);
    
    % normalised (soma-background)
    [Norm.(sess).n_bgsub,    ~]       = normalize_and_smooth( reshape(Seg.(sess).r-Seg.(sess).background,[nTm,nTr,nROI]), reshape(Seg.(sess).t,[nTm,nTr,nROI]), p );
    Norm.(sess).n_bgsub   = reshape(Norm.(sess).n_bgsub, [nTm*nTr,nROI]);

    % normalized background
    [Norm.(sess).background, Norm.(sess).background_F0] = normalize_and_smooth( reshape(Seg.(sess).background,[nTm,nTr,nROI]), reshape(Seg.(sess).t,[nTm,nTr,nROI]), p );
    Norm.(sess).background   = reshape(Norm.(sess).background, [nTm*nTr,nROI]);

    
    Norm.(sess).t   = Seg.(sess).t;
    Norm.(sess).p = p;
    % remove uncorrected frames
    Norm.(sess).n(Seg.(sess).error,:) = NaN;
    Norm.(sess).n_bgsub(Seg.(sess).error,:) = NaN;
    Norm.(sess).background(Seg.(sess).error,:) = NaN;
    MC_error.(sess) = Seg.(sess).error;
end

masks = Seg.masks; im=Seg.mean_image;
save([exp.analysed, '/masked_norm_test.mat'], 'Norm','p','masks','im','MC_error','-v7.3');

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

