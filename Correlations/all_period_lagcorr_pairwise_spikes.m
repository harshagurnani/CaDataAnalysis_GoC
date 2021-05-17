function [lags, corrs, whisk_period, varargout ] =  all_period_lagcorr_pairwise_spikes( POI, spikes, sample_time, smooth_window, ds_rate, lag_halfrange, varargin )
%% SYNTAX: [l, c] = all_period_lagcorr_pairwise_dff( locoperiod, spike (integral), sample_times,  Downsampling_rate, MaxLag [, extratime, min_POI_length, n_shuffle, do_plot )
%   Computes pairwise Correlation coefficient between spike trains for all ROIs
%   Correlation is computed after both data has been cropped to respective
%   periods, downsampled to 'ds_rate' and within a range of
%   [-lag_halfrange, lag_halfrange] (in ms).
% 
%   INPUTS:
%       - POI   [2-D array]
%           Periods Of Interest. 
%           Start and end times (ms) of periods of interest should be first
%           two columns.
%       - dff    [2-D array]
%           All dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - dff_time    [2-D array]
%           All Timestamps (ms) for dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - ds_rate   [INT]
%           Downsampling rate for computing correlations (in Hz)
%       - lag_halfrange     [INT]
%           Maximum lag(or lead) (in ms)  for cropping correlations
%
%   (optional)
%       - Varargin{1}: extratime     [INT or [INT,INT]]
%           Extratime (in ms) around period of interest.
%           Input empty if setting to default
%       - varargin{2}: min_POI_length [INT]
%           minimum size of motion period (without extratime) 
%           Input empty if setting to default
%       - varargin{3}: N_shuffle?        [INT]
%           Number of shifted (shuffle) traces used to compute null
%           distribution of correlation coefficient. Significance is
%           computed only if this value is greater than 0.
%       - varargin{4}: plot?        [BOOL]
%           Plot correlation coefficients for each period? True by default
%
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [4-D array]
%           Lag correlation values ('coeff'). Dim 1 & 2 is ROI number, Dim
%           3 is lag value (in ms). Dim 4 is period of interest.
%
%
%%%%%%%%% TO ADD!!!!
%%%%%%%%% Allow inputs to be cell arrays instead of vectors
%%%%%%%%% Allow sending cropped data
%%%%%%%%% Different modes of correlations

    %% Parse inputs
    extratime = []; min_POI_length = 0; do_plot = true; n_shuffle = 0;
    extra_args = numel(varargin);
    if extra_args > 0, extratime = varargin{1};         end
    if extra_args > 1, min_POI_length = varargin{2};    end
    if extra_args > 2, n_shuffle = varargin{3};         end
    if extra_args > 3, do_plot = varargin{3};           end
    
    if isempty(extratime),      extratime = 0;      end     % No additional padding by default
    if isempty(min_POI_length), min_POI_length = 0; end     % No minimum POI length by default
    if isempty(n_shuffle),      n_shuffle = 0;      end     % Shuffle control/significance testing is off by default

    %Corr only between +/- 1.5 s BY DEFAULT
    if isempty(lag_halfrange), lag_halfrange = 1500;    end
    
    whisk_period = [];
    %% Which whisking periods to analyse?
    for kk=1:size(POI,1)
        if POI(kk,2) - POI(kk,1) >= min_POI_length  %Minimum size of POI
            whisk_period=[whisk_period; kk];
        end
    end

    n_periods = max(size(whisk_period));    %No. of whisking periods
    n_ROIs = size(spikes,2);   %No. of ROIs

    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    if mod(size(lags,2), 2)==0
        disp('Making number of lags odd, and ensuring lag 0 included')
        lags = (-lag_halfrange-500/ds_rate:1000/ds_rate:lag_halfrange+500/ds_rate);
    end
    n_lags = floor(size(lags,2)/2);

    %Array indexing: corrs( roi1, roi2, lag, whisk_period ) 
    corrs = nan( n_ROIs,  n_ROIs, 2*n_lags+1, n_periods);
    
    tic
    dt0 = sample_time(2,1) - sample_time(1);
    smooth_bins = ceil(smooth_window/dt0);
    firingRate = spikes; for jj=1:n_ROIs, firingRate(:,jj) = smooth(spikes(:,jj), smooth_bins); end
    for roi =1:n_ROIs
        [ dff_per{roi}, dff_pertime{roi} ] = get_TSOI( firingRate , sample_time(:,roi), POI(whisk_period,:), extratime);
    end
    
    for roi1= 1:n_ROIs
    for roi2 = roi1+1:n_ROIs 
       for ww=1:max(size(whisk_period))
          try
            [c,l] = get_lagcorr_ds( dff_per{roi1}{ww}, dff_pertime{roi1}{ww}, dff_per{roi2}{ww}, dff_pertime{roi2}{ww}, ds_rate);
            %Mid = zero lag index
            mid=(size(l,2)+1)/2;
          
            %Crop between relevant lags
            if 2*n_lags+1 <= size(c,1)
                corrs(roi1, roi2,:,ww) = c( mid-n_lags:mid+n_lags);
            else 
                corrs(roi1, roi2, (n_lags+1)-(mid-1):(n_lags+1)+(mid-1), ww) = c;
            end
          catch ME
            disp(sprintf('ROI %d, ROI %d, Period %d \n %s', roi1, roi2, whisk_period(ww), ME.message))
          end
       end
    end
    end
toc

   %% Do Shuffle control? 
if n_shuffle > 0
            
    total_L = size(spikes,1);
    dt = sample_time(2,1)-sample_time(1,1);
    minshift = floor(1500/dt);   %min shift by 1.5s
    temp_corrs = nan(n_shuffle, n_ROIs, n_ROIs, n_periods);
    
    
    % interpolate before shuffling - as we are using same timeseries
    parfor shuffle_id = 1:n_shuffle
        
        tcorrs = nan( n_ROIs, n_ROIs, n_periods);
        temp_dff = nan( total_L, n_ROIs); dff_per = cell(n_ROIs, 1); dff_pertime = dff_per;
        
        % Randomly shift all traces
        for roi= 1:n_ROIs
            temp_dff(:,roi) = smooth(spikes(randperm(length(spikes)),roi), smooth_bins); %randomly distribute spikes and smooth  
            %%% Crop dff data (or second data series)
            [ dff_per{roi}, dff_pertime{roi} ] = get_TSOI( temp_dff(:,roi), sample_time(:,roi), POI(whisk_period,:), extratime);
        end
        for roi1 = 1:n_ROIs
        for roi2 = roi1+1:n_ROIs
        for ww=1:n_periods
            try
                [tcorrs(roi1, roi2, ww),~] = get_zcorr_ds( dff_per{roi1}{ww}, dff_pertime{roi1}{ww}, dff_per{roi2}{ww}, dff_pertime{roi2}{ww}, ds_rate);           
            catch
            end
        end
        end
        end
        temp_corrs( shuffle_id, :, :, : ) =tcorrs;
    end
    
    
    shuff_null_levels = nan(n_ROIs, n_ROIs, 4, n_periods );
    for roi1=1:n_ROIs
    for roi2=roi1+1:n_ROIs
    for ww=1:n_periods
        shuff_null_levels( roi1, roi2, 1, ww) = prctile(temp_corrs(:,roi1,roi2,ww),5);
        shuff_null_levels( roi1, roi2, 2, ww) = prctile(temp_corrs(:,roi1,roi2,ww),50);
        shuff_null_levels( roi1, roi2, 3, ww) = prctile(temp_corrs(:,roi1,roi2,ww),95);
        shuff_null_levels( roi1, roi2, 4, ww) = mean(temp_corrs(:,roi1,roi2,ww));
        % SIGNIFICANC = 1 if higher than 95% perctile,
        %               -1 if less than 5% prctile,
        %               0 otherwise
        corr_sig(roi1,roi2,:,ww) = 1*( corrs(roi1,roi2,:,ww) > shuff_null_levels(roi1,roi2,3,ww) ) -1*( corrs(roi1,roi2,:,ww) < shuff_null_levels( roi1,roi2,1,ww ) );          
    end
    end
    end
    if nargout > 0, varargout{1} = shuff_null_levels; end
    if nargout > 1, varargout{2} = corr_sig;   end

end
    %% Plotting
    
if do_plot
    
    
    figure()
    colormap jet
    imagesc( nanmean(corrs(:,:, n_lags+1,:),4), [-1,1] )     
    xlabel('ROI 1 '); ylabel('ROI 2');
    title('Average zero-lag correlations across all periods')  

    n_rows=min(4,ceil(sqrt(n_periods)));
    n_cols=min(4,ceil(sqrt(n_periods)));
    nsubfigs = n_rows*n_cols;
    for id=1:n_periods
        if mod(id-1,nsubfigs) == 0 
            figure();
            colormap jet
            suptitle('Zero-lag correlations in different whisking/locomotion periods')   
        end
        subplot(n_rows,n_cols, mod(id-1,9)+1)
        imagesc( corrs(:,:,n_lags+1, id) , [-1,1] )
        title(['POI ' num2str(whisk_period(id))])
        if mod(id-1,nsubfigs) == 0,                            ylabel('ROI 2'); end
        if mod(id-1,nsubfigs) == nsubfigs-1 || id==n_periods,  xlabel('ROI 1'); end
        
    end
       
end
end