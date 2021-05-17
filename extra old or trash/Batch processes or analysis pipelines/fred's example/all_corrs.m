function [lags, corrs] =  all_corrs( Whiskers, TimeAxon_whisk_pos, Axon_dFF_smooth_whisk_pos, ds_rate, lag_halfrange )
%% SYNTAX: [l, c] = all_corrs( Whisker_structure, dff_times_whiskpos, dff_values_whiskpos, Downsampling_rate, MaxLag )
%   Computes Correlation coefficient between whisker angle and dFF for all
%   axons and whisking periods. Correlatiopn is computed after both data
%   has been downsampled to 'ds_rate' and within a range of
%   [-lag_halfrange, lag_halfrange] (in ms). 
%
%   ARGUMENTS:
%       - Whiskers    [STRUCTURE}
%           Assumes the following structure (for Fred's data) with whisker
%           angle in its field 'whisking_traces_pos' and timestamps in its
%           field 'whisking_traces_pos_time'
%       - TimeAxon_whisk_pos    [2-D array]
%           Timestamps of dFF data. Row index is Whisking period, column
%           index is Axon/ROI number
%       - Axon_dFF_smooth_whisk_pos    [2-D array]
%           dFF data. Row index is Whisking period, column index is
%           Axon/ROI number
%       - ds_rate   [INT]
%           Downsampling rate for computing correlations (in Hz)
%       - lag_halfrange     [INT]
%           Maximum lag(or lead) (in ms)  for cropping correlations
%
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [3-D array]
%           Lag correlation values ('coeff'). Dim 1 is Axon/ROI number, Dim
%           2 is lag value (in ms), Dim 3 is whisking period/reindexed.
%
%


    
    whisk_period = [];
    %% Which whisking periods to analyse?
    for kk=1:max(size(Whiskers.whisking_traces_pos))
        if ~isnan(Whiskers.whisking_traces_pos{kk}) & (kk ~= 12) & (kk ~= 24)   %Some problem with time data in period 12 and 24
            whisk_period=[whisk_period; kk];
        end
    end

    %Sampling rate for computing correlation
    %ds_rate=50; %Hz

    n_periods = max(size(whisk_period));    %No. of whisking periods
    n_ROIs = size(TimeAxon_whisk_pos,2);   %No. of axons

    %CCorr only between +/- 1.5 s BY DEFAULT
    %lag_halfrange = 1000;   %ms
    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    n_lags = (size(lags,2)-1)/2;

    %Array indexong: corrs( axon_number, lag, whisk_period ) 
    corrs = nan( n_ROIs, size(lags,2), n_periods);

    for ww=1:max(size(whisk_period))
       w_prd = whisk_period(ww);

       %Read whisker position and timestamps
       whisk_time = Whiskers.whisking_traces_pos_time{w_prd};
       whisk_pos = Whiskers.whisking_traces_pos{w_prd};

       for axon = 1:n_ROIs
          %Read dFF for all axons and timestamps
          dff =  Axon_dFF_smooth_whisk_pos{w_prd, axon};
          dff_time = TimeAxon_whisk_pos{w_prd, axon};
          %Get lag correlation
          [c,l] = get_lagcorr_ds( whisk_pos, whisk_time, dff, dff_time, ds_rate);
          %Mid = zero lag index
          mid=(size(l,2)+1)/2;
          %Crop between relevant lags
          corrs(axon,:,ww) = c( mid-n_lags:mid+n_lags);       
       end


    end


    %Subplot parameters
    n_rows=floor(sqrt(n_periods));
    n_cols=ceil(n_periods/n_rows);
    n_ticks=5;  %No. of x-axis ticks
    last_tick = 2*n_lags+1;
    tick_jump = floor(last_tick/n_ticks);

    figure();
    colormap jet
    suptitle('Different whisking/locomotion periods')
    for id=1:n_periods
        subplot(n_rows,n_cols, id)
        imagesc( corrs(:,:,id) , [-1,1] )
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
    end
    
    
    n_rows=floor(sqrt(n_ROIs));
    n_cols=ceil(n_periods/n_ROIs);
    figure();
    colormap jet
    suptitle('Different ROIs')
    for id=1:n_ROIs
        subplot(n_rows,n_cols, id)
        imagesc( permute(squeeze(corrs(id,:,:)),[2,1]) , [-1,1] )       %Reshaping to Whisker_period x Lag
        set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
    end
    
end