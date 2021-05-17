function [ corrs, lags ] = get_SCC( spk_prob, times, lag, ds_rate, varargin ) 

    if iscell(times), times = cell2mat(times); end
    if iscell(spk_prob), spk_prob = cell2mat(spk_prob); end
    spk_prob = double(spk_prob);
    if isempty(varargin)
        TTot = max(times(:));
        POI = [0 TTot];
    else
        POI = varargin{1};
    end
    
    [lags, corrs, ~, ~, ~] =  pairwise_period_lagcorr( POI, spk_prob, times, ds_rate, lag, 'do_plot', false );
end