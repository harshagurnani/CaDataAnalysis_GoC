function [lags, corrs, beh_period, varargout ] =  all_period_lagcorr_pairwise_dff_new( POI, dff, dff_time, ds_rate, lag_halfrange, varargin )
%% SYNTAX: [l, c] = all_period_lagcorr_pairwise_dff( locoperiod, dff_values, dff_times,  Downsampling_rate, MaxLag [, extratime, min_POI_length, n_shuffle, do_plot, do_concat )
%   Computes pairwise Correlation coefficient between dff for all ROIs
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
%       - varargin{5}: concatenate?     [BOOL]
%           Concatenate all periods of interest - shows a warning flag if
%           asking for non-zero lags
%
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [4-D array]
%           Lag correlation values ('coeff'). Dim 1 & 2 is ROI number, Dim
%           3 is lag value (in ms). Dim 4 is period of interest.
%       - beh_period [1-D array]
%           Indices of accepted/analysed behavioural periods - based on
%           minimum length
%
%   (optional)
%       - Varargout{1}: shuff_null_levels     [4D array]
%           ROI x ROI x Dist Level x BEH Period
%           5th, 50th, 95th quantiles of the shuffled distribution as well
%           as the mean, for each pair of ROI (and each period, if not
%           concatenated).
%           Dim 1 and 2 are ROI1 and ROI2, Dim3 is which quantile, and Dim
%           4 is period number. In dim 3, index 1 is 5th percentile
%           correlation coefficient, index 2 is 50th, 3rd is 95th, and 4th
%           is the mean.
%           Eg: 
%           shuff_null_levels( r1, r2, 2, p) gives the 50th percentile
%           i.e median of the shuffle-null distribution of correlation
%           coefficients for period p, and ROIs r1 and r2. produced by
%           shuffling timeseries in period p, n_shuffle times, and
%           computing crosscorrelation for each shuffle.
%           Similarly, shuff_null_levels( r1, r2, 3, p) will be the 95th
%           percentile.
%       - varargout{2}: corr_sig [3D array]
%           Significance value of computed correlation coefficient based on
%           the shuffle-null distribution.
%           Same size as corrs array - Dim 1 and 2 is ROI 1 and 2, Dim 3 is
%           lag, Dim 4 is period - corr_sig gives significance of the
%           corresponding value in corrs. Corr sig is 1 if actuall corr
%           coeff is larger than 95th of control, -1 if it is smaller than 
%           5th of control, and 0 if not significant i.e. falls within
%           5-95 percentile of the shuffle control distribution.
%
%
%%%%%%%%% TO ADD!!!!
%%%%%%%%% Allow inputs to be cell arrays instead of vectors
%%%%%%%%% Allow sending cropped data
%%%%%%%%% Different modes of correlations

    %% Parse inputs
    extratime = []; min_POI_length = 0; do_plot = true; n_shuffle = 0; do_concat =false;
    extra_args = numel(varargin);
    if extra_args > 0, extratime = varargin{1};         end
    if extra_args > 1, min_POI_length = varargin{2};    end
    if extra_args > 2, n_shuffle = varargin{3};         end
    if extra_args > 3, do_plot = varargin{4};           end
    if extra_args > 4, do_concat = varargin{5};         end
    
    if isempty(extratime),      extratime = 0;      end     % No additional padding by default
    if isempty(min_POI_length), min_POI_length = 0; end     % No minimum POI length by default
    if isempty(n_shuffle),      n_shuffle = 0;      end     % Shuffle control/significance testing is off by default
    if isempty(do_concat), do_concat = false; end % Don't concatenate periods of interest by default
    %Corr only between +/- 1.5 s BY DEFAULT
    if isempty(lag_halfrange), lag_halfrange = 1500;    end
    
    beh_period = [];
    %% Which periods of interest to analyse?
    for kk=1:size(POI,1)
        if POI(kk,2) - POI(kk,1) >= min_POI_length  %Minimum size of POI
            beh_period=[beh_period; kk];
        end
    end

    n_periods = max(size(beh_period));  %No. of periods
    n_ROIs = size(dff,2);               %No. of ROIs

    %% Which lags
    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    if mod(size(lags,2), 2)==0
        disp('Making number of lags odd, and ensuring lag 0 included')
        lags = (-lag_halfrange-500/ds_rate:1000/ds_rate:lag_halfrange+500/ds_rate);
    end
    n_lags = floor(size(lags,2)/2);

    
    %% Interpolate/Resample first to new timepoints and then crop !
    tic
    tinterval = [max(dff_time(1,:)), min(dff_time(end,:))]; % largest shared time interval
    [tmp, ~] = resample_ds( dff(:,1), dff_time(:,1), ds_rate, tinterval ); n1 = length(tmp);
    dff_new = nan(n1, n_ROIs); dff_newtime = dff_new;
    for roi =1:n_ROIs
        [dff_new(:,roi), dff_newtime(:,roi)] = resample_ds( dff(:,roi), dff_time(:,roi), ds_rate, tinterval );              % Resample
%         dff_new(:,roi) = inpaint_nans( dff_new(:,roi),4);
        [ dff_per{roi}, ~ ] = get_TSOI( dff_new(:,roi), dff_newtime(:,roi), POI(beh_period,:), extratime);                  % Crop in periods of interest
    end
    
    if do_concat
        if lag_halfrange ~= 0, error('Ask for non-zero lag correlations! Cannot concatenate trials'); end
        for roi=1:n_ROIs
        % Get concatenated traces
            dff_per{roi}{1} = cell2mat(dff_per{roi}');
        end
        n_periods = 1;
    end
    
    %Array indexing: corrs( roi1, roi2, lag, beh_period ) 
    corrs = nan( n_ROIs,  n_ROIs, 2*n_lags+1, n_periods);
    
    
    for roi1 = 1:n_ROIs
    for roi2 = roi1+1:n_ROIs 
       for ww=1:n_periods
          A1 = squeeze(dff_per{roi1}{ww}); A2 = squeeze(dff_per{roi2}{ww});
          try
            [c,l] = xcorr(A1-mean(A1), A2-mean(A2),'coeff');    %lag-correlations
            l = l*1000/ds_rate;   %ms
            
            %Mid = zero lag index
            mid=(size(l,2)+1)/2;
          
            %Crop between relevant lags
            if 2*n_lags+1 <= size(c,1)
                corrs(roi1, roi2,:,ww) = c( mid-n_lags:mid+n_lags);
            else 
                corrs(roi1, roi2, (n_lags+1)-(mid-1):(n_lags+1)+(mid-1), ww) = c;
            end
          catch ME
            disp(sprintf('ROI %d, ROI %d, Period %d \n %s', roi1, roi2, beh_period(ww), ME.message))
          end
       end
       corrs(roi2, roi1,:,:) = corrs(roi1, roi2,:,:);
    end
    end
toc

%-----------------------------------------------------------------------------
   %% Do Shuffle control? 
   tic
if n_shuffle > 0
            
    total_L = size(dff_new,1);
    dt = dff_newtime(2,1)-dff_newtime(1,1);
    minshift = floor(1500/dt);   %min shift by 1.5s
    temp_corrs = nan(n_shuffle, n_ROIs, n_ROIs, n_periods);
    
    parfor shuffle_id = 1:n_shuffle
        % NOTE: shuffle first then crop to all periods of interest
        
        tcorrs = nan( n_ROIs, n_ROIs, n_periods);
        temp_dff = nan( total_L, n_ROIs); dff_per = cell(n_ROIs, 1); dff_pertime = dff_per;
        
        % Circularly shift all traces by random, independent amounts
        for roi= 1:n_ROIs
            shift = ceil(unifrnd(minshift, total_L-minshift));
            temp_dff(:,roi) =circshift(dff_new(:,roi), shift);  
            %%% Crop dff data (or second data series)
            [ dff_per{roi}, ~ ] = get_TSOI( temp_dff(:,roi), dff_newtime(:,roi), POI(beh_period,:), extratime);
        end
        
        if do_concat
            for roi=1:n_ROIs
            % Get concatenated traces
            dff_per{roi}{1} = cell2mat(dff_per{roi}');
            end
        end
        
        for roi1 = 1:n_ROIs
        for roi2 = roi1+1:n_ROIs
            for ww=1:n_periods
                try
                A1 = squeeze(dff_per{roi1}{ww}); A2 = squeeze(dff_per{roi2}{ww});
                cc = corrcoef(A1, A2)
                tcorrs(roi1, roi2, ww) = cc(2); 
                catch
                end
            end
        end
        end
        temp_corrs( shuffle_id, :, :, : ) =tcorrs;
    end
    
    
    shuff_null_levels = nan(n_ROIs, n_ROIs, 4, n_periods ); % saves 5th, 50th,95th percentile and mean of shuffle distribution instead of complete distribution
    corr_sig = nan(size(corrs));
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
    shuff_null_levels( roi2, roi1, :, :) = shuff_null_levels( roi1, roi2, :, :);
    end
    end
    if nargout > 0, varargout{1} = shuff_null_levels; end
    if nargout > 1, varargout{2} = corr_sig;   end

end
    toc
    
%-----------------------------------------------------------------------------   
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
        title(['POI ' num2str(beh_period(id))])
        if mod(id-1,nsubfigs) == 0,                            ylabel('ROI 2'); end
        if mod(id-1,nsubfigs) == nsubfigs-1 || id==n_periods,  xlabel('ROI 1'); end
        
    end
       
end
end