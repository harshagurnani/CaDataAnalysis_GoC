%% Prepare data

% Load data - load p,t,r,n. 
if ~exist( 'r', 'var') || ~exist('p', 'var') 
    [p,r,t,~] = load_POI();
end

nroi = size(r, 3); tm = size(r, 1); trl = size(r, 2);

% Getting only the good parameters
if exist('g', 'var')
    if isstruct(g)
        [all_traces_raw, time_traces_raw] = avg_grouped_roi( g, r, t, p, 'donorm', false );
        nsoma = length(g.soma);
        nog = length(g.other_groups);
        x = cell( 1, nsoma+nog); xt = x;
        for ss = 1:nsoma, x{ss} = reshape( all_traces_raw{ss}(:,:,1), [tm*trl, 1] ); xt{ss} = reshape( time_traces_raw{ss}(:,:,1), [tm*trl, 1] ); end
        for gg = 1:nog, x{nsoma+gg} = reshape( all_traces_raw{end-nog+gg}(:,:,1), [tm*trl, 1] ); xt{nsoma+gg} = reshape( time_traces_raw{end-nog+gg}(:,:,1), [tm*trl, 1] ); end
        find_gid = false;
        xt = cell2mat(xt);
    else 
        find_gid=true;
    end
else
    find_gid = true;
end

% Find good ROI with high maximal dff values
if find_gid
    if ~exist('gid','var')
    zscored = (r - repmat( mean(mean(r,1),2), [tm, trl, 1] ) )./(repmat(reshape( std( reshape(r, [tm*trl, nroi]) ) ,[1 1 nroi] ),[tm,trl,1]));
    pid = prctile( reshape(zscored, [tm*trl, nroi]), 95, 1 );
    gid = find( pid >= 0.5*prctile(pid,80) );       % Atleast half as bright as the top 20% of ROIs
    end
    x = reshape( r(:,:,gid), [tm*trl, length(gid)] ); xt = reshape( t(:,:,gid), [tm*trl, length(gid)] );
    x = mat2cell( x, tm*trl, ones(1, length(gid)));    
end

nROI_select = numel(x);
% for jj=1:nROI_select, x{jj} = x{jj}/prctile(x{jj},10);  end         %Use 10 percentile as baseline rather than the mean?

%% Spike Inference

% Set the drift parameter
driftrate_perc = 0.01;

% Set the estimation algorithm
estimate = 'map'; % Spike probabilities. 
%Use MAP for a single spike train - the map estimate. or 'Samples' for returning multiple spike trains?


% Use rise time for GCamp6f
% rise_time = 15;%ms

% Run algorithm from Deneux et al, 2016
run_method = 'batch';   % same variance for observation noise
[spk, fit, drift, parest]= deneux_spike_inference( x,   1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method );
% [spk fit drift parest]= deneux_spike_inference( x,
% 1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method, 'rise_time', risetime );