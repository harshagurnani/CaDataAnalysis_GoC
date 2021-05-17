function [ data, infer ] = benchmark_multiple_roi( p, r, t, g, varargin)
%% Benchmark spike inference without ground truth data - Infer spike trains and compute cross correlations for multiple pixels on same cell/group
% Usage: 1) [data, infer] = benchmark_multiple_roi( p, r, t, g )
%               Pools traces and infers spikes for each group. Sets all
%               options to default
%        2) [data, infer] = benchmark_multiple_roi( p, r, t, g , opts )
%               Uses additional options based on fields in opts
%        3) [data, infer] = benchmark_multiple_roi( p, r, t, g , 'batch' )
%                        OR benchmark_multiple_roi( p, r, t, g , 'individual' )
%               Sets inference method to 'batch' (or 'individual') for both
%               pooled and individual data
%        4) [data, infer] = benchmark_multiple_roi( p, r, t, g , false )
%               Does not do inference/plotting/correlations. Only fields in
%               data and infer are cell numbers and pooled traces
%        5) [data, infer] = benchmark_multiple_roi( p, r, t, g , false(or opts.do_infer=false), data_old, infer_old )
%               Does not perform inference but uses spike trains in
%               infer_old for plotting and correlations
%
%   Inputs:
%   ---------------------------
%       p:  [Struct]    Standard dat loading params. Only needs to have field
%                       acquisition rate (in Hz) to calculate dt
%       r:  [3d matrix] Raw fluorescence ( Time x Trial x Pixel/Patch/Mask )
%       t:  [3d matrix] Imaging times    ( Time x Trial x Pixel/Patch/Mask )
%       g:  [Struct]    Fields: 'soma', 'dend', 'other_groups'. Each field
%                       stores a cell array of different ROI. a cell
%                       contains a vector of Pixel ids corresponding to
%                       that ROI. Eg. for 2 soma: g.soma = 1x2 cell array.
%                       g.soma{1} = [1 5] where pixel=1,5 belong to soma 1
%
%   (Optional inputs):
%       opts: [Struct] - Options for inference (see usage above)
%               opts.do_infer     [Logical]  - Do spike inference?
%               opts.group_method [Char]     - 'batch' or 'individual' for Noise estimation for pooled data
%               opts.ind_method   [Char]     - 'batch' or 'individual' for Noise estimation for individual pixel data
%               opts.do_plot      [Logical]  - Plot inferred spike trains?
%               opts.do_corr      [Logical]  - Do cross-correlation for groups of spiketrains?
%
%       infer_old:  Used to compute correlations or plotting if no spike
%                   inference performed
%       data_old:   Returned as data if no spike inference performed
%
%   Outputs:
%   ----------------------------
%       data:   [Struct]    data.group (Pooled) and data.ind (individual)
%                           data (with subfields 'raw' and 'times')
%                           data.ind.group_id has group numbers for all
%                           individual pixels
%       infer:  [Struct]    infer.group (Pooled) and infer.ind (individual)
%                           with subfields: 'spikes', 'drift', 'params'
%                           (estimated spike inference parameters.
%                           infer.corr stores cross-correlations for groups
%                           with subfields 'smooth_bin', 'lags',
%                           'crosscorr'
%
%   Harsha G. 
%   Nov 2017
% ------------------------------------------------------------------------%



%% Parsing inputs
opts.do_infer = true;
opts.group_method = 'individual';
opts.ind_method = 'batch';
opts.do_plot = true;    opts.do_corr = true;
data_given = false;

opts_given = false;
if nargin>4
    if isstruct( varargin{1} )
        tmp = varargin{1};
        if isfield( tmp, 'group_method'), opts.group_method = tmp.group_method; end
        if isfield( tmp, 'ind_method'),   opts.ind_method   = tmp.ind_method;   end
        if isfield( tmp, 'do_infer'),     opts.do_infer     = tmp.do_infer;     end
        if isfield( tmp, 'do_plot'),      opts.do_plot      = tmp.do_plot;      end
        if isfield( tmp, 'do_corr'),      opts.do_corr      = tmp.do_corr;      end
        opts_given = true;
    elseif islogical( varargin{1} )
        opts.do_infer = varargin{1};
        opts_given = true;
    elseif ischar( varargin{1} )
        opts.group_method = varargin{1};
        opts.ind_method   = varargin{1};
        opts_given = true;
    end
end

if nargin>5
    if opts_given, data_id = 2; else, data_id = 1; end
    data = varargin{data_id}; infer = varargin{data_id+1};
    data_given = true;
end

if ~data_given && ~opts.do_infer
    opts.do_plot = false; opts.do_corr= false;
end

%% Pooling data in groups
tm = size(r, 1); trl = size(r, 2);

% Creating the averaged traces
[all_traces_raw, time_traces_raw] = avg_grouped_roi( g, r, t, p, 'donorm', false );
nsoma = length(g.soma);
nog   = length(g.other_groups);
pooled_f = cell( 1, nsoma+nog); pooled_t = pooled_f;

for ss = 1:nsoma, pooled_f{ss} = reshape( all_traces_raw{ss}(:,:,1), [tm*trl, 1] ); pooled_t{ss} = reshape( time_traces_raw{ss}(:,:,1), [tm*trl, 1] ); end
for gg = 1:nog, pooled_f{nsoma+gg} = reshape( all_traces_raw{end-nog+gg}(:,:,1), [tm*trl, 1] ); pooled_t{nsoma+gg} = reshape( time_traces_raw{end-nog+gg}(:,:,1), [tm*trl, 1] ); end
pooled_t = cell2mat(pooled_t);


% Reordering traces based on groups
new_id = []; cell_number = [];
for ss = 1:nsoma, nid = length(g.soma{ss}(:));          new_id = [new_id;  g.soma{ss}(:) ];           cell_number = [cell_number,         ss*ones(1,nid)]; end
for gg = 1:nog
    group = g.other_groups{gg}; 
    if iscell(group), group = group{1}; end 
    nid = length(group(:));  new_id = [new_id;  group(:) ];   cell_number = [cell_number, (nsoma+gg)*ones(1,nid)]; 
end
r = r(:,:, new_id); t = t(:,:,new_id);
nroi = length(new_id);

% Find good individual ROI with high maximal dff values belonging to the
% cells above

zscored = (r - repmat( mean(mean(r,1),2), [tm, trl, 1] ) )./(repmat(reshape( std( reshape(r, [tm*trl, nroi]) ) ,[1 1 nroi] ),[tm,trl,1]));
pid = prctile( reshape(zscored, [tm*trl, nroi]), 95, 1 );
include_id = find( pid >= 0.5*prctile(pid,80) );       % Atleast half as bright as the top 20% of ROIs

ind_f = reshape( r(:,:,include_id), [tm*trl, length(include_id)] ); ind_t = reshape( t(:,:,include_id), [tm*trl, length(include_id)] );     cell_number = cell_number(include_id);
ind_f = mat2cell( ind_f, tm*trl, ones(1, length(include_id)));    


%% Spike Inference for Averaged Groups 
if opts.do_infer
 
    % Set the drift parameter
    driftrate_perc = 0.01;

    % Set the estimation algorithm
    estimate = 'map'; % Spike probabilities. 
    %Use MAP for a single spike train - the map estimate. or 'Samples' for returning multiple spike trains?


    % Use rise time for GCamp6f
    % rise_time = 15;%ms

    % Run algorithm from Deneux et al, 2016
    run_method = opts.group_method;   % different observation noise because different number of ROI may have been averaged.

    sprintf('Running spike inference on %d averaged traces... ZZZ ...', numel(pooled_f))
    [spk_pooled, fit_pooled, drift_pooled, parest_pooled]= deneux_spike_inference( pooled_f,   1/p.acquisition_rate, 'est_method', estimate, 'driftrate', driftrate_perc, 'run_method', run_method );
    % [spk fit drift parest]= deneux_spike_inference( x,
    % 1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method, 'rise_time', risetime );

    
%% Spike Inference for Individual ROIs 

    % Set the drift parameter
    driftrate_perc = 0.01;

    % Set the estimation algorithm
    estimate = 'map'; % Spike probabilities. 
    %Use MAP for a single spike train - the map estimate. or 'Samples' for returning multiple spike trains?


    % Use rise time for GCamp6f
    % rise_time = 15;%ms

    % Run algorithm from Deneux et al, 2016
    run_method = ind.group_method;   % different observation noise because different number of ROI may have been averaged.

    sprintf('Running spike inference on %d individual scanned point ROI... ZZZ ...', numel(ind_f) )
    [spk_ind, fit_ind, drift_ind, parest_ind]= deneux_spike_inference( ind_f,   1/p.acquisition_rate, 'est_method', estimate, 'driftrate', driftrate_perc, 'run_method', run_method );
    % [spk fit drift parest]= deneux_spike_inference( x,
    % 1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method, 'rise_time', risetime );

else
    
    
    
end


if data_given &&  ~opts.do_infer
    spk_pooled =infer.group.spikes;
    spk_ind    = infer.ind.spikes;
    pooled_t   = data.group.time;
    ind_t      = data.ind.time;
end

%% Plotting spike trains
if opts.do_plot
    nGroups = length(pooled_f);
    nROI    = numel(ind_f);
    c=nROI+nGroups*3+3;

    figure;hold on;
    for gg = 1:nGroups
       c=c-3;
       % Plot spike train from averaged trace
       spike_raster = spk_pooled{gg}; if iscell(spike_raste), spike_raster = spike_raster{1}; end
       spike_raster(spike_raster == 0) = NaN;
       plot( pooled_t(:,gg), spike_raster+c, 'LineWidth',2 );
       c=c-1;
       % Plot spike trains from individual imaged poimt ROI
       group_roi = find(cell_number==gg);
       for jj=1:length(group_roi)
           indx = group_roi(jj); spike_raster = spk_ind{indx}; 
           spike_raster(spike_raster == 0) = NaN;
           plot( ind_t(:,indx), spike_raster+c, 'LineWidth',1 );
           c=c-1;
       end

    end
end

%% Calculating total cross-correlation within groups
if opts.do_corr
    tmax = max(ind_t(:));
    smooth_bin = 100;    %ms
    dt = 1000/p.acquisition_rate;  %ms
    [ lags, lag_crosscorr ] = calc_SCC( spk_pooled, spk_ind, pooled_t, ind_t, cell_number, smooth_bin, dt, tmax );
end

%% Preparing returned args

data.group.raw = pooled_f;  data.group.time = pooled_t;
data.ind.raw   = ind_f;     data.ind.time   = ind_t;
data.ind.group_id = cell_number;

if opts.do_infer
    infer.group.spikes = spk_pooled;   infer.group.drift = drift_pooled;    infer.group.params = parest_pooled;
    infer.ind.spikes   = spk_ind;      infer.ind.drift   = drift_ind;       infer.ind.params   = parest_ind;
end

if opts.do_corr
    infer.corr.lags = lags;     infer.corr.crosscorr = lag_crosscorr;   infer.corr.smooth_bin = smooth_bin;
end

end

function [ lags, lag_crosscorr ] = calc_SCC( spk_pooled, spk_ind, pooled_t, ind_t, groupID, smooth_bin, dt, tmax )
% to calculate pairwise CC between spike (FR) estimates of points on same
% region-of-interest
  
  smooth_factor = 2*(floor(smooth_bin*0.5/dt))+1;
  nGroups = numel(spk_pooled);
  for gg = 1:nGroups
      indx = ( groupID == gg );
      avg_trace = spk_pooled{gg}; if iscell(avg_trace), avg_trace = avg_trace{1}; end
      x  = [avg_trace, cell2mat(spk_ind(indx))];  %concatenate group-and-ind traces
      xt = [pooled_t(:,gg), ind_t(:,indx)];
      x=double(x);
      for nr = 1:size(x,2)
          x(:,nr) = smooth(x(:,nr), smooth_factor );
      end
      [lags, group_SCC{gg}, ~, ~, ~ ] =  pairwise_period_lagcorr( [0 tmax] , x, xt, 1000/smooth_bin, 500, 'do_plot', false );
      if iscell(group_SCC{gg}), group_SCC{gg}=group_SCC{gg}{1}; end
  end
    nlags = size(group_SCC{1},3); nR = nGroups+numel(spk_ind); lag_crosscorr = nan(nR, nR, nlags);
    c=1;
    for gg=1:nGroups
          nrows= sum(groupID==gg);
          lag_crosscorr(c:c+nrows,c:c+nrows,:) = group_SCC{gg};
          c=c+nrows+1;
%         lag_crosscorr(:,:,ll) = blkdiag(group_SCC{:});
    end
end