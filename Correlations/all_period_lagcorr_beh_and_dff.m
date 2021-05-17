function [lags, corrs, beh_period, varargout] =      all_period_lagcorr_beh_and_dff( POI, viddata, vidtimes, dff, dff_time, ds_rate, lag_halfrange, varargin )
%% SYNTAX: [l, c] = all_period_lagcorr( locoperiod, speed, speedtime, dff_values_whiskpos, dff_times_whiskpos,  Downsampling_rate, MaxLag [, extratime, min_POI_length, n_shuffle, do_plot )
%   Computes Correlation coefficient between viddata and dff for all ROIs
%   and whisking/locomotion periods. For eg. between whisker angle and
%   fluorescence. 
%   Correlation is computed after both data has been cropped to respective
%   periods, downsampled to 'ds_rate' > Correlations are kepts for lags
%   within a range of [-lag_halfrange, lag_halfrange] (in ms).
%   Optionally, significance can be tested by computing a 'null
%   distribution' of correlation coefficients for each ROI and period. Corr
%   is considered significant if outside the [5,95] range of the
%   corresponding null distribution.
%
%   NOTE: Normalisation method for xcorr is 'coeff'(lag-zero autocorr is 1), and
%   all traces are mean-subtracted before calculating Corr coeffs
%
%   INPUTS:
%       - POI   [2-D array]
%           Periods Of Interest. 
%           Start and end times (ms) of periods of interest should be first
%           two columns.
%       - vidtimes    [1-D array]
%           Timestamps (ms) for 'viddata'
%       - viddata    [1-D array]
%           Motion index or speed or some behavioural data etc. Naming is
%           only representative. (Is a single timeseries that is not
%           different across ROI )
%       - dff    [2-D array]
%           All dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - dff_time    [2-D array]
%           All Timestamps (ms) for dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - ds_rate   [INT]
%           Downsampling rate for computing correlations (in Hz)
%       - lag_halfrange     [INT]
%           Maximum lag(or lead) (in ms)  for cropping(keeping) correlation
%           matrix(vector)
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
%
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [3-D array]
%           Lag correlation values ('coeff'). Dim 1 is ROI number, Dim
%           2 is lag value (in ms), Dim 3 is period of interest
%       - beh_period [1-D array]
%           List of period of interest that were analysed and plotted.
%
%   (optional)
%       - varargout{1}: shuff_null_levels [3-D array]
%           Different percentile levels in the null distribution computed
%           for each ROI and whisk period. Dim 1 is ROI number, Dim 2 has 4
%           columns corresponding to 4 levels, Dim 3 is period of interest.
%           The 4 levels are percentile 5, percentile 50(median),
%           percentile 95 and the mean.
%       - varargout{2}: corrs_sig         [3-D Boolean array]
%           Significance value for each correlation coefficient in corrs. 
%           Value is 1 if it is higher than 95th percentile of null dist, 
%           -1 if it is lower than 5th percentile of null dist, and
%           0 otherwise.
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
    
    beh_period = [];
    %% Which behaviour periods to analyse?
    for kk=1:size(POI,1)
        if POI(kk,2) - POI(kk,1) >= min_POI_length  %Minimum size of POI
            beh_period=[beh_period; kk];
        end
    end
    n_periods = max(size(beh_period));
    
    %Sampling rate for computing correlation
    %ds_rate=50; %Hz

    n_periods = max(size(beh_period));    %No. of whisking periods
    n_ROIs = size(dff,2);   %No. of ROIs

    
    %lag_halfrange = 1000;   %ms
    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    if mod(size(lags,2), 2)==0
        disp('Making number of lags odd, and ensuring lag 0 included')
        lags = (-lag_halfrange-500/ds_rate:1000/ds_rate:lag_halfrange+500/ds_rate);
    end
    n_lags = floor(size(lags,2)/2);
    

    %Array indexing: corrs( axon_number, lag, beh_period ) 
    corrs = nan( n_ROIs, 2*n_lags+1, n_periods);


    tic
    %%% Crop video data (or first data series)
    [ whisk_mi, whisk_time ] = get_TSOI(viddata, vidtimes, POI(beh_period,:), extratime);
    
    for roi= 1:n_ROIs
%         dff(:,roi) = inpaint_nans( dff(:,roi),4);
        %%% Crop dff data (or second data series)
        [ dff_per, dff_pertime ] = get_TSOI( dff(:,roi), dff_time(:,roi), POI(beh_period,:), extratime);
        for ww=1:n_periods
           
          %Get lag correlation
          try
            [c,l] = get_lagcorr_ds( whisk_mi{ww}, whisk_time{ww}, dff_per{ww}, dff_pertime{ww}, ds_rate);
            %Mid = zero lag index
            mid=(size(l,2)+1)/2;
          
            %Crop between relevant lags
            if 2*n_lags+1 <= size(c,1)
                corrs(roi,:,ww) = c( mid-n_lags:mid+n_lags);       
            else 
                corrs(roi, (n_lags+1)-(mid-1):(n_lags+1)+(mid-1), ww) = c;
            end
          catch ME
            disp(sprintf('ROI %d, Period %d \n %s', roi, beh_period(ww), ME.message))
          end
       end


    end
    toc

    
    %% Do shuffle control??
if n_shuffle > 0
            
    totalL = size(dff,1);
    dt = dff_time(2,1)-dff_time(1,1);
    minshift = floor(1500/dt);   %min shift by 1.5s

    temp_corrs = nan(n_shuffle, n_ROIs, n_periods);
    
    parfor shuffle_id = 1:n_shuffle
    tcorrs = nan( n_ROIs, n_periods);
    for roi= 1:n_ROIs
        shift = ceil(unifrnd(minshift, totalL-minshift));
        temp_dff =circshift(dff, shift, 1);  
        %%% Crop dff data (or second data series)
        [ dff_per, dff_pertime ] = get_TSOI( temp_dff(:,roi), dff_time(:,roi), POI(beh_period,:), extratime);
        for ww=1:n_periods
          %Get correlation coeff for mismatched traces
          try
            [tcorrs(roi, ww),~] = get_zcorr_ds( whisk_mi{ww}, whisk_time{ww}, dff_per{ww}, dff_pertime{ww}, ds_rate);           
          catch ME
            disp(sprintf('ROI %d, Period %d \n %s', roi, beh_period(ww), ME.message))
          end
       end


    end
    temp_corrs( shuffle_id, :, : ) =tcorrs;
    end
    toc

    shuff_null_levels = nan(n_ROIs, 4, n_periods );
    for roi=1:n_ROIs
        for ww=1:n_periods
            shuff_null_levels( roi, 1, ww) = prctile(temp_corrs(:,roi,ww),2.5);
            shuff_null_levels( roi, 2, ww) = prctile(temp_corrs(:,roi,ww),50);
            shuff_null_levels( roi, 3, ww) = prctile(temp_corrs(:,roi,ww),97.5);
            shuff_null_levels( roi, 4, ww) = mean(temp_corrs(:,roi,ww));
            %1 if higher than 95% perctile, -1 if less than 5% prctile, 0 otherwise
            corr_sig(roi,:,ww) = 1*( corrs(roi,:,ww) > shuff_null_levels(roi,3,ww) ) -1*( corrs(roi,:,ww) < shuff_null_levels( roi,1,ww ) );          
        end
    end
    
    if nargout > 0, varargout{1} = shuff_null_levels; end
    if nargout > 1, varargout{2} = corr_sig;   end
    
end

    
    %% Plotting
if do_plot

    n_ticks=7;  %No. of x-axis ticks
    last_tick = 2*n_lags+1;
    tick_jump = floor(last_tick/n_ticks);    
    
    % Plot 1 - Lag correlations for each ROI. Each sublot is a different
    % period of interest
    
    %Subplot parameters
    n_rows=min(4,ceil(sqrt(n_periods)));
    n_cols=min(5,ceil(sqrt(n_periods)));
    nsubfigs = n_rows*n_cols;
    for id=1:n_periods
        if mod(id-1,nsubfigs) == 0 
            figure();   %new figure;
            colormap jet
            suptitle('Different periods')   
        end
        subplot(n_rows,n_cols, mod(id-1,20)+1)
        imagesc( corrs(:,:,id) , [-1,1] )
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis#
        title(['POI ' num2str(beh_period(id))])
        if mod(id-1,nsubfigs) == 0,                            ylabel('ROI number'); end
        if mod(id-1,nsubfigs) == nsubfigs-1 || id==n_periods,  xlabel('Lag (ms)'); end
        
    end
    

    % Plot 2 - Lag correlation in each period. Each subplot is a different
    % ROI
    
    n_rows=min(5,ceil(sqrt(n_ROIs)));
    n_cols=n_rows;
    nsubfigs = n_rows^2;
    for id=1:n_ROIs
        if mod(id-1,nsubfigs) == 0 
            figure();
            colormap jet
            suptitle('Different ROIs')  
        end
        subplot(n_rows,n_cols, mod(id-1,25)+1 )
        imagesc( permute(squeeze(corrs(id,:,:)),[2,1]) , [-1,1] )       %Reshaping to Whisker_period x Lag
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
        if mod(id-1,nsubfigs) == 0,                           ylabel('Period of interst'); end
        if mod(id-1,nsubfigs) == nsubfigs-1|| id==n_periods,  xlabel('Lag (ms)'); end
        
    end
    
end


end
