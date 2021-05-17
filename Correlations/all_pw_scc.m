function [ corrs ] = all_pw_scc( spk_prob, times, lag, ds_rate, beh_periods ) 
%% Spike count correlations - over entire recording window and during different behavioural periods

% Pairwise correlations - Load lp, mip, l_mi_p, no_l_mi_p etc etc into
% beh_periods

if isempty(beh_periods)
    beh_periods = struct;
end


nRates = length(ds_rate);

% Correlation at rates mentioned in ds_rate

%% over entire traces
for jj = 1:nRates
    corrs(jj).binsize_ms = 1000/ds_rate(jj); corrs(jj).lags=[];
    [ corrs(jj).full, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj) );
end

%% During running epochs
if isfield( beh_periods, 'lp' )
for jj = 1:nRates
    [ corrs(jj).running.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.lp );
    corrs(jj).running.mean = squeeze(nanmean(corrs(jj).running.individual,4));
end
end
 
%% During whisking
if isfield( beh_periods, 'mip' )
for jj = 1:nRates
    [ corrs(jj).whisking.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.mip );
    corrs(jj).whisking.mean = squeeze(nanmean(corrs(jj).whisking.individual,4));
end
end

%% During whisking and no running
if isfield( beh_periods, 'no_l_mip' )
for jj = 1:nRates
    [ corrs(jj).whisking_noloco.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.no_l_mip );
    corrs(jj).whisking_noloco.mean = squeeze(nanmean(corrs(jj).whisking_noloco.individual,4));
end
end

%% During air puffs
if isfield( beh_periods, 'puff' )
for jj = 1:nRates
    [ corrs(jj).air_puff.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.puff );
    corrs(jj).air_puff.mean = squeeze(nanmean(corrs(jj).air_puff.individual,4));
end
end

%% During self whisking and no running
if isfield( beh_periods, 'no_l_nopuff_mip' )
for jj = 1:nRates
    [ corrs(jj).self_whisking_noloco.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.no_l_nopuff_mip );
    corrs(jj).self_whisking_noloco.mean = squeeze(nanmean(corrs(jj).self_whisking_noloco.individual,4));
end
end


%% During no whisking and no running - QUIET epochs
if isfield( beh_periods, 'quiet' )
for jj = 1:nRates
    [ corrs(jj).quiet.individual, corrs(jj).lags ] = get_SCC( spk_prob, times, lag, ds_rate(jj), beh_periods.quiet );
    corrs(jj).quiet.mean = squeeze(nanmean(corrs(jj).quiet.individual,4));
end
end

end