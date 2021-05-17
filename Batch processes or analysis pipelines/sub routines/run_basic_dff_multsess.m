function [ Basic, Seg, Norm ] = run_basic_dff_multsess( exp, Basic )
    
    Basic = basic_dff_params( Basic );
    
    if strcmpi(Basic.ROI_type, 'points')
    
    else
        %% Patches and Planes
        %  do MC and segmentation on patches and planes
        
        do_ref  = true;
        % setup ref for movement correction - if already specified
        if ~isempty( exp.Use.Ref )
            load( exp.Use.Ref, 'MC' );
            ref = MC.ref;   ref_roi = MC.ref_ROI; do_ref = false; clear MC
            MC = struct( 'ref', ref, 'ref_ROI', ref_ROI' );
            clear ref ref_roi
        end

        %% For all sessions - Load Data + Motion Correction
        for sn=1:exp.nSess

        sess            = exp.session_types{exp.use_sessions(sn)};

        % Load data
        [p,r,t,~] = load_POI('source',exp.(sess),     'ROI_type', Basic.ROI_type, 'ROIs', Basic.ROIs, 'trials', Basic.trials,...
                             'avg_patch',Basic.patch_avg,  'normalize', false );
        [nL, nP, nTm, nTr, nROI ] = size(r);

        if ~Basic.continuous_time && nTr>1                                  % all trial times are identical
            t(:,:,:,2:nTr, :) = repmat( t(:,:,:,1,:), [1,1,1,nTr-1,1]);
        end

        if Basic.MC                                                         % do Movement correction
            MC.max_disp = Basic.MC_max_disp;
            run_MC_pipeline;
        end

        if numel(size(r))==5                                                % Reshaoe to 4 D : [x, y, time, roi]
            r = reshape(r, [nL, nP, nTm*nTr, nROI]);        t= reshape( t, [nL, nP, nTm*nTr, nROI] );
        end

        % Save data
        Basic.(sess) = struct( 'r',r,'t',t,'params',p );
        save( [exp.analysed, '\raw_', sess, '.mat'], 'Basic', 'MC', 'exp', 'p', '-v7.3' );

        % save mean image and pixels present
        if ~exist('im','var')
            im = cell(nROI,1); for jj=1:nROI, im{jj} = nan(nL,nP,exp.nSess); end
            keeppix = im;
        end
        for roi=1:nROI
            im{roi}(1:nL, 1:nP, sn)      = nanmean(r(:,:,:,roi),3);           % mean image
            keeppix{roi}(1:nL, 1:nP, sn) = ~isnan(mean(r(:,:,:,roi),3));      % pixels present in all frames
        end


        disp(['.....................Saved Imaging Data for ', sess, ' session after motion correction'])

        % show result
        if isempty(Basic.MC_show), MC_show = strcmpi(input('Do you want to test MC quality? y/n:    ','s'), 'y'); end
        if Basic.MC_show 
                nframes = min(50, 0.05*nTm*nTr);	% show random N frames
                for tm=sort(unidrnd(nTm*nTr, 1,nframes))
                for roi=1:nROI, subplot(nr,nr,roi);  imagesc(r(:,:,tm,roi),[50 2000]);   end
                suptitle(['Frame ',num2str(tm)]);    pause(0.002)
                end
        end
        clear r t MC_show 
        
        end
        for sn=1:exp.nSess
            sess            = exp.session_types{exp.use_sessions(sn)};
            Basic.(sess)    = [];
        end
        %% ----------------------------------------------------------------------

        %% Masks /Segmentation
        clear MC

        for roi=1:nROI
        im{roi}      = nanmean(im{roi}, 3);                  % mean image across all sessions
        keeppix{roi} = (sum(keeppix{roi},3) == exp.nSess) ;      % pixels present in all sessions
        end
        Seg.mean_image = im;


        get_masks_multsess;

        %% Apply masks
        for sn = 1:exp.nSess
            sess            = exp.session_types{exp.use_sessions(sn)};
            load( [exp.analysed, '\raw_', sess, '.mat'], 'Basic', 'MC' );
            apply_masks_multsess;
            clear Basic MC
        end
        save( [exp.analysed, '\masked_raw.mat'], 'Seg', 'exp', 'keeppix', '-v7.3' );
        disp('.....................Saved Masked Data')
        clear errorid keeppix

    end

    
    
    %% Normalisation
    % in paint NaNs only for normalization process

    for sn = 1:exp.nSess
        sess            = exp.session_types{exp.use_sessions(sn)};
        [~, nROI ] = size(Seg.(sess).r);
        for roi=1:nROI
            Seg.(sess).r(:,roi) = inpaint_nans( Seg.(sess).r(:,roi),4);  
        end

        % use defaults for now
        load( [exp.analysed, '\raw_', sess, '.mat'], 'p' );
        load( [exp.analysed, '\Params_basic.mat' ]);
        p.normalize = true;        
        p.smooth_causal = Basic.smooth_causal;
        p.smooth_scale = Basic.smooth_scale;    %ms %default
        p.base_percentile = Basic.base_percentile;
        p.baseline_pertrial = Basic.baseline_pertrial;
        Norm.(sess).p = p;        
        nTm = p.nTimepoints; nTr = size(Seg.(sess).t,1)/nTm;
        nROI = size(Seg.(sess).t, 2 );
        
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

        % remove uncorrected frames
        Norm.(sess).n(Seg.(sess).error,:) = NaN;
        Norm.(sess).n_bgsub(Seg.(sess).error,:) = NaN;
        Norm.(sess).background(Seg.(sess).error,:) = NaN;
        MC_error.(sess) = Seg.(sess).error;
        
        
    end

    masks = Seg.masks; im=Seg.mean_image;
    save([exp.analysed, '\masked_norm_test.mat'], 'Norm','p','masks','im','MC_error','-v7.3');

    clear MC_error masks im

    disp('.....................Saved Normalised Traces')

end
