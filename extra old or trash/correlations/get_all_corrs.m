whisk_period = [];
%% Which whisking periods to analyse?
for kk=1:max(size(Whiskers.whisking_traces_pos))
    if ~isnan(Whiskers.whisking_traces_pos{kk}) & (kk ~= 12) & (kk ~= 24)   %Some problem with time data in period 12 and 24
        whisk_period=[whisk_period; kk];
    end
end

%Sampling rate for computing correlation
ds_rate=50; %Hz

n_periods = max(size(whisk_period));    %No. of whisking periods
n_axons = size(TimeAxon_whisk_pos,2);   %No. of axons

%CCorr only between +/- 1.5 s BY DEFAULT
lag_halfrange = 1000;   %ms
lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
n_lags = (size(lags,2)-1)/2;

%Array indexong: corrs( axon_number, lag, whisk_period ) 
corrs = nan( n_axons, size(lags,2), n_periods);

for ww=1:max(size(whisk_period))
   w_prd = whisk_period(ww);
   
   %Read whisker position and timestamps
   whisk_time = Whiskers.whisking_traces_pos_time{w_prd};
   whisk_pos = Whiskers.whisking_traces_pos{w_prd};
   
   for axon = 1:n_axons
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
n_rows=4;
n_cols=5;
n_ticks=5;  %No. of x-axis ticks
last_tick = 2*n_lags+1;
tick_jump = floor(last_tick/n_ticks);

figure();
colormap jet
suptitle('Different whisking periods')
for id=1:n_periods
    subplot(n_rows,n_cols, id)
    imagesc( corrs(:,:,id) , [-1,1] )
    set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
end

figure();
colormap jet
suptitle('Different axons')
for id=1:n_axons
    subplot(n_rows,n_cols, id)
    imagesc( permute(squeeze(corrs(id,:,:)),[2,1]) , [-1,1] )       %Reshaping to Whisker_period x Lag
    set(gca, 'Xtick', 1:tick_jump:last_tick, 'XTickLabel',lags(1:tick_jump:last_tick))  %Re-label x-axis
end