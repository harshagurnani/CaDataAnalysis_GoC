function [lags, corrs, whisk_period] =  all_period_lagcorr( POI, viddata, vidtimes, dff, dff_time, ds_rate, lag_halfrange, varargin )
%% SYNTAX: [l, c] = all_period_lagcorr( locoperiod, speed, speedtime, dff_values_whiskpos, dff_times_whiskpos,  Downsampling_rate, MaxLag [, extratime )
%   Computes Correlation coefficient between viddata and dff for all ROIs
%   and whisking/locomotion periods. For eg. between whisker angle and
%   fluorescence. 
%   Correlation is computed after both data has been cropped to respective
%   periods, downsampled to 'ds_rate' and within a range of
%   [-lag_halfrange, lag_halfrange] (in ms).
% 
%   ARGUMENTS:
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
%           Maximum lag(or lead) (in ms)  for cropping correlations
%
%   (optional)
%       - Varargin{1}: extratime     [INT or [INT,INT]]
%           Extratime (in ms) around period of interest.
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [3-D array]
%           Lag correlation values ('coeff'). Dim 1 is ROI number, Dim
%           2 is lag value (in ms), Dim 3 is period of interest
%       - whisk_period [1-D array]
%           List of period of interest that were analysed and plotted.
%
%
%%%%%%%%% TO ADD!!!!
%%%%%%%%% Allow inputs to be cell arrays instead of vectors
%%%%%%%%% Allow sending cropped data
%%%%%%%%% Different modes of correlations

    if nargin == 8 
        extratime = varargin{1};
    else
        extratime = 300;
    end
    whisk_period = [];
    %% Which whisking periods to analyse?
    for kk=1:size(POI,1)
        if POI(kk,2) - POI(kk,1) >= 2000-extratime(1) && POI(kk,1)>=0  %At least 50 ms long
            whisk_period=[whisk_period; kk];
        end
    end

    %Sampling rate for computing correlation
    %ds_rate=50; %Hz

    n_periods = max(size(whisk_period));    %No. of whisking periods
    n_ROIs = size(dff,2);   %No. of ROIs

    %Corr only between +/- 1.5 s BY DEFAULT
    %lag_halfrange = 1000;   %ms
    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    if mod(size(lags,2), 2)==0
        disp('Making number of lags odd, and ensuring lag 0 included')
        lags = (-lag_halfrange-500/ds_rate:1000/ds_rate:lag_halfrange+500/ds_rate);
    end
    n_lags = floor(size(lags,2)/2);
    

    %Array indexing: corrs( axon_number, lag, whisk_period ) 
    corrs = nan( n_ROIs, 2*n_lags+1, n_periods);


    tic
    %%% Crop video data (or first data series)
    [ whisk_mi, whisk_time ] = get_TSOI(viddata, vidtimes, POI(whisk_period,:), extratime);
    
    for roi= 1:n_ROIs
        
        %%% Crop dff data (or second data series)
        [ dff_per, dff_pertime ] = get_TSOI( dff(:,roi), dff_time(:,roi), POI(whisk_period,:), extratime);
        for ww=1:max(size(whisk_period))
            
            %%%% To add: what if data length < lag - don't pad with zeros!!
            %%%% compute less lags!
%            data_halfrange = POI(whisk_period(ww),2)-POI(whisk_period(ww),1)/2 + min(extratime);   %Max data range
%            max_data_lag = data_halfrange*ds_rate/1000;
%            max_lag = min(n_lags, max_data_lag);
        
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
            disp(sprintf('ROI %d, Period %d \n %s', roi, whisk_period(ww), ME.message))
          end
       end


    end
    toc

    %Subplot parameters
    n_rows=4;
    n_cols=5;
    n_ticks=7;  %No. of x-axis ticks
    last_tick = 2*n_lags+1;
    tick_jump = floor(last_tick/n_ticks);

    
    for id=1:n_periods
        if mod(id-1,20) == 0 
            figure();
            colormap jet
            suptitle('Different whisking/locomotion periods')   
        end
        subplot(n_rows,n_cols, mod(id-1,20)+1)
        imagesc( corrs(:,:,id) , [-1,1] )
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
    end
    
    
    n_rows=5;
    n_cols=5;
    for id=1:n_ROIs
        if mod(id-1,25) == 0 
            figure();
            colormap jet
            suptitle('Different ROIs')  
        end
        subplot(n_rows,n_cols, mod(id-1,25)+1 )
        imagesc( permute(squeeze(corrs(id,:,:)),[2,1]) , [-1,1] )       %Reshaping to Whisker_period x Lag
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
    end
    
end
