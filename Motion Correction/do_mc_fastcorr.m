function [ ref, save_offsets, varargout ] = do_mc_fastcorr( patches, varargin )
%% Register multiple stacks(patch timeseries) simultaneously using fast cross-correlation maximisation
%  Register each individual stack (corresponding to different ROI) to a
%  reference by finding offset at each timepoint that maximises
%  cross-correlation of shifted image with the reference. 
%  Each patch can have a different offset or a single offset can be applied
%  to all patches at the same timepoint. Multiple reference and offset
%  selection options detailed below.
%   ----------------------------------
%           USAGE: 
%   ---------------------------------
%       2 to 6 outputs:
%       - output 1 and 2 are always reference (ref) and offsets. 
%       - If 'get_corr' is true, output 3 is peak crosscorrelation for each
%       frame&roi.
%       - The successive optional outputs are registered stacks, expanded
%       reference, and params structure.
%
%    a.     [ref, offsets ] = do_mc_fastcorr( patches, 'Name', Value...) ...
%               returns reference and offsets
%                       OR
%   If 'get_corr' is false (default) :
%   ---------------
%    b.     [A, B, new_reg, new_ref ] = do_mc_fastcorr( patches, 'Name', Value...)
%               also returns the registered stack (expanded frame), and its
%               correspondingly expanded reference.
%                       OR
%    c.     [A, B, C, D, params ] = do_mc_fastcorr( patches, 'Name', Value...)
%               also returns the offset and reference selection parameters as structure object
%                       OR
%   If 'get_corr' is true :
%   ---------------
%    d.     [A, B, peakCC ] = do_mc_fastcorr( patches, 'Name', Value...)
%               also returns peak crosscorrelation
%                       OR
%    e.     [A, B, peakCC, new_reg, ... ] = do_mc_fastcorr( patches, 'Name', Value...)
%               returns peak crosscorrelation and other optional outputs
%
%   ----------------------------------
%           INPUTS:
%   ----------------------------------
%       - patches   [3 or 4 or 5D array] 
%                   Intensities (raw or normalized) of frame at each
%                   timepoint for all ROI.
%                   Dim 1 and 2 are frame x-y, Dim 3 is time. A single ROI timeseries can be
%                   supplied as 3D or 4D array. If array is 5D, dim 4 is
%                   read as trial and dim 5 as ROI number, and array is
%                   reshaped into 4D(trials concatenated).
%
%   (optional) Supplied as Name-Value Pairs
%           To specify reference as MATRIX, or SLICE NUMBERS TO AVERAGE:
%           -----------
%       - 'ref'     To set up the reference for registration.
%                   1) If value is a vector, it is read as frame numbers
%                   (dim 3) to be averaged to create reference.
%                   2) If value is a matrix (x-y-roi) or a cell array(
%                   ref{roi}), then corresponding frame is used as
%                   reference. x-y size must be same size as patches data.
%
%           TO create reference by finding PERIODS OF NO MOTION (SPEED=0
%           for MIN_PERIOD) and use 30 FRAMES from REF_START time of the
%           first such period:
%       - 'speed'       Speed/motion timeseries (2-column array). Column 1 is
%                       time, column 2 is speed. Can be used to determine frame
%                       numbers to use for creating reference image. Must be
%                       used in conjuction with arguments 'times'. 
%       - 'times'       [Vector] Timestamp for each frame (same length as total
%                       number of frames for a single ROI i.e TimexTrial)
%       - 'min_period'  [NUM] Time in ms - min length of no motion period.
%                       5000 by default.
%       - 'ref_start'   [NUM] Time in ms - time from start of no motion
%                       period  to use for reference. 
%                       1500 by default.
%
%       If none of the options are specified, or no silent period is found,
%       the FIRST 30 FRAMES are used to create reference
%
%           Offset Selection:
%           -----------
%       - 'opt_offset'  [BOOL] If true(default), a single optimal offset is
%                       determined for all ROI at each timepoint. 
%       - 'opt_method'  [STR] Method to find optimal offset across all ROI
%                       Can be one of 'weighted COM' (default), 'COM', or
%                       'valid ROI COM'
%       - 'all_offsets' [BOOL] If true, individual offsets are also
%                       returned in the params structure. False by default
%
%
%       - 'use_offsets' [2D or 3D array] X & Y offsets for each timepoint
%                       (and roi)
%
%
%           Additional:
%           -----------
%       - 'get_corr'    [BOOL] If true, returns peak crosscorrelation
%                       values for each timepoint and roi as 2D array
%
%       - 'do_crop'     [BOOL] If true, crops registered image into same
%                       size as original patch data, rather than returning
%                       expanded frame
%
%       - 'verbose'     [BOOL] If true(default), execution times for
%                       different sections of the code will be printed.
%
%   ----------------------------------
%           OUTPUTS:
%   ----------------------------------
%       - ref       Cell array (ROIx1). Reference image used to register
%                   each ROI's stack 
%       - offsets   [2D(t x 2) or 3D(t x roi x 2) array] X & Y offsets that
%                   maximize cross correlation between reference and
%                   functional image at each timepoint. Different offsets
%                   for each ROI if 'opt_offset' set to false.
%
%   (optional) If you want the offsets to be applied to return registered
%   stacks
%       - peakCC    2D array (t-roi)
%                   Peak cross correlation
%       - new_reg   4D array (x-y-t-roi)
%                   Registered stack with expanded frames (rather than
%                   cropping).
%       - new_ref   Cell array with the same reference as in ref but with
%                   expanded frames (same frame size as new_reg)
%       - p         Structure with offset and reference options
%
%
%   Author:  Harsha G. 
%   Created: April 2018
%------------------------------------------------------------------------------------------------------------------------------------------------%    
    
    %% Parse arguments

    nLines = size(patches,1);   nPixels = size(patches,2);
    nr = nLines; nc= nPixels;
    % convert to 4-D array (lines x pixels per line x all time (trials are concatenated) x ROI 
    if length(size(patches))==5,    patches = reshape(patches, [nLines, nPixels, size(patches,3)*size(patches,4), size(patches,5)] ); end
    nTBins = size(patches, 3);  nROI = size(patches, 4);

    [p, ref, speed, t ] = parse_arguments( varargin{:} );
    
    % If offsets are given, simply apply them.
    if isfield(p.offset, 'given'),  do_mc = false;
    else,   do_mc = true;    end
        
    if do_mc
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    %%      FFT of Mean subtracted ROI
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    
    patch_fft = cell(nROI,1);
    % patch_fft is a cell array of cell arrays (convert to mat later?)
    % patch_fft{roi}{frame} = 2-D array of fourier coefficients
    tic;
    parfor roi=1:nROI
        patch_fft{roi} = run_fft( squeeze(patches(:,:,:,roi)), 2*nLines-1, 2*nPixels-1, nTBins );
    end
    t1=toc;
    if p.verbose, fprintf( ' Fourier transformed %d patches with %d frames in %f seconds\n', nROI, nTBins, t1); end
    
    
    %---------------------------------------------------------------------------------------------------------------------------------------------%    
    %%      Setting up reference for registration
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    
    % ref/ref_fft are cell arrays 
    % ref{roi} is the reference for each ROI to register to
    % ref_fft is its Fourier decomposition

    tic
    if      p.ref.given
        % convert to cell array if not
        if ~iscell( ref ),  ref = squeeze(mat2cell( ref, 2*nLines-1, 2*nPixels-1, ones(nROI,1)));   end
    elseif  p.ref.use_slicenum
        % create refs with slice numbers
        ref = arrayfun( @(roi) nanmean( patches(:,:,p.ref.ref_slices, roi),3), 1:nROI, 'UniformOutput', false );        
    elseif  p.ref.use_movedata
        % find period of zero movement(locomotion or whiking) and use 30 frames in its middle by default
        speed_dt = (speed(2,1)-speed(1,1));  bin5s = min( floor(p.min_period/speed_dt), length(speed) ); 
        sumspd = arrayfun( @(jj) sum(abs(speed(jj+1:jj+bin5s,2) )), 0:length(speed)-bin5s );    
        noloco = find(sumspd == 0, 1);
        if isempty(noloco)
        	p.ref.ref_slices = 1:30;    warning('No 5s silent periods found, using first 30 frames to create registration reference')
            ref = arrayfun( @(roi) nanmean( patches(:,:,p.ref.ref_slices, roi),3), 1:nROI, 'UniformOutput', false );      
        else
            reftime = p.ref_start + speed_dt*noloco;   %choose ref_start time (1.5s by default) into the noloco period for reference
            ref_slices = find(t>reftime,1); p.ref.ref_slices = (ref_slices:ref_slices+30);
            ref = arrayfun( @(roi) nanmean( patches(:,:,p.ref.ref_slices, roi),3), 1:nROI, 'UniformOutput', false );       
        end
        clear speed_dt sumspd noloco reftime
    else
        p.ref.ref_slices = 1:30;
        ref = arrayfun( @(roi) nanmean( patches(:,:,p.ref.ref_slices, roi),3), 1:nROI, 'UniformOutput', false );        
    end
    ref_fft = arrayfun( @(roi) fft2(ref{roi} - mean2(ref{roi}), 2*nLines-1, 2*nPixels-1), 1:nROI, 'UniformOutput', false );

    
    t2=toc;
    if p.verbose, if ~p.ref.given, fprintf( 'Set up registration reference for %d ROI in %f seconds\n', nROI, t2);   end;   end
    
    
    %---------------------------------------------------------------------------------------------------------------------------------------------%    
    %%      Calculate cross-correlation
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    offsets2 = nan( nTBins, nROI, 2 );
    switch p.method
        case 'fft'
            nX = 2*nLines-1; nY = 2*nPixels-1;
            tic;
            [CC, maxCC ] = deal(cell(nROI,1));
            peakCC = nan( nTBins, nROI );
            for roi=1:nROI
                CC{roi} = nan( 2*nLines-1, 2*nPixels-1, nTBins );
                [autoCC, maxCC{roi}] = deal(nan(nTBins,1));                                    
                refCC = ifft2(ref_fft{roi}.*conj(ref_fft{roi}));   refCC = refCC(1);                     % ref auto correlation
                for jj=1:nTBins
                    CC{roi}(:,:,jj) = fftshift( ifft2(ref_fft{roi}.*conj(patch_fft{roi}{jj})) );         % total cross correlation
                    maxCC{roi}(jj) = max(max(CC{roi}(:,:,jj)));                                          % max cross correlation(unnormalised)
                    tmp= ifft2(patch_fft{roi}{jj}.*conj(patch_fft{roi}{jj})) ;    autoCC(jj) = tmp(1);   % auto correlation
                end

                % Normalise by autoCC of frame and reference
                peakCC(:,roi) = maxCC{roi}./sqrt(refCC*autoCC);                                 % max cross correlation(Normalised) - to use for combining patch offsets
            end

            t3=toc;
            if p.verbose, fprintf( ' Calculated normalised cross correlations in %f seconds\n', t3);    end
            clear patch_fft
            
        case 'dft'
            nX = nLines; nY = nPixels;
            parfor roi=1:nROI
                patch_fft{roi} = run_fft( squeeze(patches(:,:,:,roi)), nLines, nPixels, nTBins );
            end
            ref_fft = arrayfun( @(roi) fft2(ref{roi} - mean2(ref{roi}), nLines, nPixels), 1:nROI, 'UniformOutput', false );

            
            tic;
            [CC, maxCC ] = deal(cell(nROI,1));
           
            peakCC = nan( nTBins, nROI );
            for roi=1:nROI
                CC{roi} = zeros(  nLines, nPixels, nTBins );
                
                
                [autoCC, maxCC{roi}] = deal(nan(nTBins,1));                                    
                refCC = ifft2(ref_fft{roi}.*conj(ref_fft{roi}));   refCC = refCC(1);                     % ref auto correlation
                for jj=1:nTBins
                    [output] = dftregistration(patch_fft{roi}{jj},ref_fft{roi},p.usfac);
                    xshift = -output(4);                    yshift = -output(3);
                    offsets2(jj, roi,:) = [xshift, yshift];
                    CC{roi}(:,:,jj) = fftshift( ifft2(ref_fft{roi}.*conj(patch_fft{roi}{jj})) );                      % total cross correlation
                    maxCC{roi}(jj) = max(max(CC{roi}(:,:,jj)));                                          % max cross correlation(unnormalised)
%                     autoCC(jj) = 1;   % auto correlation
                end

                % Normalise by autoCC of frame and reference
                peakCC(:,roi) = maxCC{roi};                                 % max cross correlation(Normalised) 
            end

            t3=toc;
            if p.verbose, fprintf( ' Calculated normalised cross correlations in %f seconds\n', t3);    end
            clear patch_fft
    end
    
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    %%      Find optimal offsets
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    tic;
    
    % find offsets that maximize CC for each patch
    if p.method =='fft'
        for roi = 1:nROI
        s = nan(nTBins,1,2);
        for jj=1:nTBins
            [s(jj,1,1), s(jj,1,2)] = ind2sub( [nX,nY], find(CC{roi}(:,:,jj)== maxCC{roi}(jj),1) );
        end
        offsets2(:,roi,:) = s;
        end
    
        offsets2 = offsets2 - repmat(reshape([nLines nPixels],[1 1 2]), [nTBins, nROI, 1]); %centred 
    end
    
    if p.offset.optimize
    % Single optimal offset across patches
        offsets = nan(nTBins, 2);
        if strcmp(p.offset.method, 'weighted COM' )
            % Peak CC is used as weight
            if p.method=="fft"
                parfor jj=1:nTBins
                offsets( jj, : ) = peakCC(jj,:) * squeeze(offsets2(jj,:,:)) / sum(peakCC(jj,:)) ;
                end
            elseif p.method=="dft"
                parfor jj=1:nTBins
                newVec = [peakCC(jj,:) * exp(1i*2*pi*offsets2(jj,:,1)/nr)', ...
                          peakCC(jj,:) * exp(1i*2*pi*squeeze(offsets2(jj,:,2))/nc)'] / sum(peakCC(jj,:));
                offsets(jj,:) = [nr,nc].*[angle(newVec)]/(2*pi)
                end 
            end
        elseif strcmp(p.offset.method, 'COM' )
            % COM of offsets from all roi
            parfor jj=1:nTBins
                offsets( jj, : ) = nanmean(squeeze(offsets2(jj,:,:)));
            end
            
        elseif strcmp(p.offset.method, 'valid ROI COM' )
            % COM of all roi that have a valid/good peak
             
        end
        if p.return_all_offsets
            p.offsets.all = offsets2;
        end    
    else
        offsets = offsets2;
        
    end
    
 
    t4=toc;
    if p.verbose, fprintf( ' Optimised offsets in %f seconds\n', t4); end
    
    else
       offsets = p.offset.given;
       new_ref = [];
    end
    if p.return_corr,   p.corr = peakCC;    end
    
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    %%      Apply offsets (optional)
    %---------------------------------------------------------------------------------------------------------------------------------------------%
    apply_corr = false; rega = 1; 
    if ~p.return_corr
        if nargout>2, apply_corr = true;    end
    elseif nargout == 3
       varargout{1} = peakCC;
    elseif nargout>3
        varargout{1} = peakCC;
        apply_corr = true;   rega = 2;
    end
    save_offsets = offsets; % Return raw non-integral offsets
    
    if apply_corr
        if p.method =='fft'
            % correct offsets befor applying!
            % 1) blank offsets that are more than maxo (if set)
            if ~isempty(p.offset.maxo),  offsets( abs(offsets) > p.offset.maxo )=NaN;  end
            % 2) interpolate offsets that are NaN
            tmp = reshape(offsets, [size(offsets,1), size(offsets,2)*size(offsets,3) ]);
            for jj=1:size(tmp,2), tmp(:,jj) =inpaint_nans(tmp(:,jj)); end
            offsets = reshape(tmp, size(offsets,1), size(offsets,2), size(offsets,3));
            if size(offsets,3)==1, offsets=squeeze(offsets); end
            % 3) Re-clamp offsets
            if ~isempty(p.offset.maxo)
                offsets( offsets >  p.offset.maxo ) =  p.offset.maxo;  
                offsets( offsets < -p.offset.maxo ) = -p.offset.maxo;  
            end
            % offsets are integral (no subpixel registration)
            offsets = real(round(offsets)); 

            if length(size(offsets))==2
                %same offset to all

                crop = [ max( 0, max(-offsets(:,1))), max(0, max(offsets(:,1))) ;
                         max( 0, max(-offsets(:,2))), max(0, max(offsets(:,2))) ];

                     
                if strcmpi(p.newreg_default , 'nans')
                    new_reg = nan( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2),   nTBins,  nROI );
                else
                    new_reg = zeros( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2),   nTBins,  nROI ); % nan -> zeros, to allow successive registration
                end
                for jj=1:nTBins
                    new_reg( 1+crop(1,1)+offsets(jj,1):nLines+crop(1,1)+offsets(jj,1), 1+crop(2,1)+offsets(jj,2):nPixels+crop(2,1)+offsets(jj,2), jj, : )= patches(:,:,jj,:) ;
                end

                if do_mc || p.ref.given
                nROI_ref = numel(ref);
                new_ref = cell(nROI_ref,1);
                for roi=1:nROI_ref
                    new_ref{roi} = nan( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2) );
                    new_ref{roi}(1+crop(1,1):end-crop(1,2), 1+crop(2,1):end-crop(2,2)) = ref{roi};
                end
                end



            else
                %different offsets to each patch
                crop = [ max( 0, max(max(-offsets(:,:,1)))), max(0, max(max(offsets(:,:,1)))) ;
                         max( 0, max(max(-offsets(:,:,2)))), max(0, max(max(offsets(:,:,2)))) ];

                if strcmpi(p.newreg_default , 'nans')
                    new_reg = nan( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2),   nTBins,  nROI );
                else
                    new_reg = zeros( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2),   nTBins,  nROI );
                end
                for roi=1:nROI
                for jj=1:nTBins
                    new_reg( 1+crop(1,1)+offsets(jj,roi,1):nLines+crop(1,1)+offsets(jj,roi,1), 1+crop(2,1)+offsets(jj,roi,2):nPixels+crop(2,1)+offsets(jj,roi,2), jj, : )= patches(:,:,jj,:) ;
                end
                end

                if do_mc || p.ref.given
                    nROI_ref = numel(ref);
                    new_ref = cell(nROI_ref,1);
                    for roi=1:nROI_ref
                        new_ref{roi} = nan( nLines+crop(1,1)+crop(1,2),  nPixels+crop(2,1)+crop(2,2) );
                        new_ref{roi}(1+crop(1,1):end-crop(1,2), 1+crop(2,1):end-crop(2,2)) = ref{roi};
                    end
                    if p.crop
                         new_ref = arrayfun( @(roi) new_ref{roi}(crop(1,1)+1:nLines+crop(1,1), crop(2,1)+1:nPixels+crop(2,1), :, : ) ,1:nROI_ref, 'UniformOutput', false);
                    end
                end


            end
            if p.crop
               new_reg = new_reg(crop(1,1)+1:nLines+crop(1,1), crop(2,1)+1:nPixels+crop(2,1), :, : );

            end
        
        elseif p.method =='dft'
            if size(offsets,3)==1
                % apply same offsets to all
                Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
                [Nc,Nr] = meshgrid(Nc,Nr);
                new_reg = nan( nr, nc, nTBins, nROI);
                for roi=1:nROI
                for jj=1:nTBins
                    new_reg(:,:,jj,roi) = ifft2(patch_fft{roi}{jj}.*exp(1i*2*pi*(-offsets(jj,1)*Nr/nr-offsets(jj,2)*Nc/nc)));
                end
                end
            else
                % different offsets
                Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
                [Nc,Nr] = meshgrid(Nc,Nr);
                new_reg = nan( nr, nc, nTBins, nROI);
                for roi=1:nROI
                for jj=1:nTBins
                    new_reg(:,:,jj,roi) = ifft2(patch_fft{roi}{jj}.*exp(1i*2*pi*(-offsets(jj,roi,1)*Nr/nr-offsets(jj,roi,2)*Nc/nc)));
                end
                end
            end
        end
        
        if nargout>=rega+2
            varargout{rega} = new_reg;
        end
        if nargout>rega+2
            varargout{rega+1} = new_ref;
        end
        if  nargout>2+(rega+1)
            varargout{rega+2} = p;
        end
    
    end

    
end




%-----------------------------------------------------------------------%
% HELPER

function fw = run_fft( data, nx, ny, nt )
    fw = arrayfun( @(frame) fft2(data(:,:,frame) - mean2(data(:,:,frame)), nx, ny ), 1:nt, ...
                                   'UniformOutput', false );  %Cell array of all patches in (Fourier)freq domain
end

function [p, ref, speed, t ] = parse_arguments( varargin )
    
    [p.ref.given, p.ref.use_slicenum, p.ref.use_movedata]  = deal(false);
    [ref, speed, t ] = deal([]);
    p.min_period = 5000;        % in ms - min size of "silent/no movement period"
    p.ref_start = 1500;         % in ms - what frame to start the reference at, after the start of the silent period?
    p.offset.optimize = true;
    p.offset.method = 'weighted COM';
    p.return_all_offsets = false;
    p.verbose = false;
    p.crop = false;
    p.return_corr = false;
    p.offset.maxo = [];
    p.method = 'fft';   %or dft
    p.newreg_default = 'nans';
    
    usfac = 1;  %no upsampling for now - qqqqqq
    p.usfac = usfac; 
    
    c=1;
     while c<=nargin
         if strcmp( varargin{c}, 'ref' )
              if isvector(varargin{c+1}) && ~iscell(varargin{c+1})
                  p.ref.use_slicenum = true;
                  p.ref.ref_slices = varargin{c+1};
              else
                ref = varargin{c+1};
                p.ref.given = true;
              end
              c = c+2;
         elseif strcmp( varargin{c}, 'speed' )
              p.ref.use_movedata = true;
              speed = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'times' )
              t = varargin{c+1};
              c = c+2;     
         elseif strcmp( varargin{c}, 'min_period' )
              p.min_period = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'ref_start' )
              p.ref_start = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'opt_offset' )
              p.offset.optimize = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'opt_method' )
              p.offset.optimize = true;
              p.offset.method = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'all_offsets' )
              p.return_all_offsets = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'verbose' )
              p.verbose = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'use_offsets' )
              p.offset.given = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'do_crop' )
              p.crop = varargin{c+1};
              c = c+2;
         
         elseif strcmp( varargin{c}, 'get_corr' )
              p.return_corr = varargin{c+1};
              c = c+2;
         
         elseif strcmp( varargin{c}, 'max_disp' )
              p.offset.maxo = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'method' )
              p.method = varargin{c+1};
              c = c+2;
         
         elseif strcmp( varargin{c}, 'usfac' )
              p.usfac = varargin{c+1};
              c = c+2;
         elseif strcmp( varargin{c}, 'newreg_default' )
              p.newreg_default = varargin{c+1};
              c = c+2;
         else 
              error('No argument by the name %s', varargin{c} )
         end
       
    end
    
end