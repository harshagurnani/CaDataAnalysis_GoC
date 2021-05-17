%% Apply masks
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    suffix          = sess;
%     load( [exp.analysed, '\raw_', suffix, '.mat'], 'Basic', 'MC' );

    [ Seg.(sess).r,  Seg.(sess).t ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.masks);           % avg somatic fluorescence
    [ Seg.(sess).background, ~ ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.background_masks);   % avg dark baclground fluorescence
    
    % remove uncorrected frames
    errorid = find(MC_error.(sess));
    Seg.(sess).error = MC_error.(sess);
    if max(diff(errorid))<8, Seg.(sess).error(errorid(1):errorid(end)) = true; end
    Seg.(sess).r( Seg.(sess).error,: ) = NaN;
    
    
end

for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    suffix          = sess;
    for roi=1:nROI
        Seg.(sess).r(:,roi) = inpaint_nans( Seg.(sess).r(:,roi),4);  
    end

    % use defaults for now (100 ms smoothing with moving average)
%     load( [exp.analysed, '\raw_', suffix, '.mat'], 'p' );
    
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
%     MC_error.(sess) = Seg.(sess).error;
end
