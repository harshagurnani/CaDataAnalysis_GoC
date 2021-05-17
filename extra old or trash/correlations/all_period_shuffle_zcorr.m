function [ corrs, corr_sig, corr_rho, shuff_corr, whisk_period] =  all_period_shuffle_zcorr( POI, viddata, vidtimes, dff, dff_time, ds_rate, n_shuffle, varargin )
%% SYNTAX: [l, c] = all_period_lagcorr( locoperiod, speed, speedtime, dff_values_whiskpos, dff_times_whiskpos,  Downsampling_rate, num_shuffles [, extratime, min_shuffletime )
%   Computes Correlation coefficient between viddata and dff for all ROIs
%   and whisking/locomotion periods. For eg. between whisker angle and
%   fluorescence. 
%   Correlation is computed after both data has been cropped to respective
%   periods, downsampled to 'ds_rate'. The rho of the linear correlation,
%   as well as the 5,50,95 percentile and mean of correlations with
%   circularly shuffled behaviour data.
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
%       - n_shuffle     [INT]
%           Number of shuffles
%
%   (optional)
%       - Varargin{1}: extratime     [INT or [INT,INT]]
%           Extratime (in ms) around period of interest.
%   OUTPUT:

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

   

    %Array indexing: corrs( axon_number, lag, whisk_period ) 
    corrs = nan( n_ROIs, n_periods);
    corr_sig = zeros( n_ROIs, n_periods);
    corr_rho = corrs;
    shuff_corr = nan( n_ROIs, n_periods, 4 );

    tic
    %%% Crop video data (or first data series)
    [ whisk_mi, whisk_time ] = get_TSOI(viddata, vidtimes, POI(whisk_period,:), extratime);
    
    %% Correlation coefficient
    for roi= 1:n_ROIs
        
        %%% Crop dff data (or second data series)
        [ dff_per, dff_pertime ] = get_TSOI( dff(:,roi), dff_time(:,roi), POI(whisk_period,:), extratime);
        for ww=1:max(size(whisk_period))
          try
            [corrs( roi, ww),corr_rho(roi,ww)] = get_zcorr_ds( whisk_mi{ww}, whisk_time{ww}, dff_per{ww}, dff_pertime{ww}, ds_rate);
          catch ME
            disp(sprintf('ROI %d, Period %d \n %s', roi, whisk_period(ww), ME.message))
          end
       end
    end
    
    totalL = size(dff,1);
    minshift = 2000/30;
    temp_corrs = nan(n_shuffle, n_ROIs, n_periods);
    
    parfor shuffle_id = 1:n_shuffle
    tcorrs = nan( n_ROIs, n_periods);
    for roi= 1:n_ROIs
        shift = ceil(unifrnd(minshift, totalL-minshift));
        temp_dff =circshift(dff, shift, 1);  
        %%% Crop dff data (or second data series)
        [ dff_per, dff_pertime ] = get_TSOI( temp_dff(:,roi), dff_time(:,roi), POI(whisk_period,:), extratime);
        for ww=1:max(size(whisk_period))
            
        
          %Get lag correlation
          try
            [tcorrs(roi, ww),~] = get_zcorr_ds( whisk_mi{ww}, whisk_time{ww}, dff_per{ww}, dff_pertime{ww}, ds_rate);
           
          catch ME
            disp(sprintf('ROI %d, Period %d \n %s', roi, whisk_period(ww), ME.message))
          end
       end


    end
    temp_corrs( shuffle_id, :, : ) =tcorrs;
    end
    toc

    for roi=1:n_ROIs
        for ww=1:max(size(whisk_period))
            shuff_corr( roi, ww, 1) = prctile(temp_corrs(:,roi,ww),5);
            shuff_corr( roi, ww, 2) = prctile(temp_corrs(:,roi,ww),50);
            shuff_corr( roi, ww, 3) = prctile(temp_corrs(:,roi,ww),95);
            shuff_corr( roi, ww, 4) = mean(temp_corrs(:,roi,ww));
            if corrs( roi, ww ) > shuff_corr( roi, ww, 3)
                corr_sig(roi,ww) = 1;
            elseif corrs(roi,ww) < shuff_corr( roi,ww,1)
                corr_sig(roi,ww) = -1;
            end
        end
    end
    
end
