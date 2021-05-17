function [ data, infer ] = benchmark_multiple_roi_commonnoise( p, r, t, g)


tm = size(r, 1); trl = size(r, 2);

% Creating the averaged traces
[all_traces_raw, time_traces_raw] = avg_grouped_roi( g, r, t, p, 'donorm', false );
nsoma = length(g.soma);
nog = length(g.other_groups);
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

% Set the drift parameter
driftrate_perc = 0.01;

% Set the estimation algorithm
estimate = 'map'; % Spike probabilities. 
%Use MAP for a single spike train - the map estimate. or 'Samples' for returning multiple spike trains?


% Use rise time for GCamp6f
% rise_time = 15;%ms

% Run algorithm from Deneux et al, 2016
run_method = 'batch';   % different observation noise because different number of ROI may have been averaged.

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
run_method = 'batch';   % different observation noise because different number of ROI may have been averaged.

sprintf('Running spike inference on %d individual scanned point ROI... ZZZ ...', numel(ind_f) )
[spk_ind, fit_ind, drift_ind, parest_ind]= deneux_spike_inference( ind_f,   1/p.acquisition_rate, 'est_method', estimate, 'driftrate', driftrate_perc, 'run_method', run_method );
% [spk fit drift parest]= deneux_spike_inference( x,
% 1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method, 'rise_time', risetime );



%% Plotting spike trains
nGroups = length(pooled_f);
nROI    = numel(ind_f);
toplot_traces = true;

c=nROI+nGroups*3+3;
if toplot_traces
    figure;hold on;
   for gg = 1:nGroups
       c=c-3;
       % Plot spike train from averaged trace
%        spike_raster = spk_pooled{gg}{1}; spike_raster(spike_raster == 0) = NaN;
       spike_raster = spk_pooled{gg}; spike_raster(spike_raster == 0) = NaN;
       plot( pooled_t(:,gg), spike_raster+c, 'LineWidth',2 );
       c=c-1;
       % Plot spike trains from individual imaged poimt ROI
       group_roi = find(cell_number==gg);
       for jj=1:length(group_roi)
           indx = group_roi(jj); spike_raster = spk_ind{indx}; spike_raster(spike_raster == 0) = NaN;
%            size(spike_raster);
           plot( ind_t(:,indx), spike_raster+c, 'LineWidth',1 );
           c=c-1;
       end
   end
end


%% Calculating total cross-correlation within groups
tmax = max(ind_t(:));
smooth_bin = 40;    %ms
dt = 1000/p.acquisition_rate;  %ms
[ lags, lag_crosscorr ] = calc_SCC( spk_pooled, spk_ind, pooled_t, ind_t, cell_number, smooth_bin, dt, tmax );
 
%% Preparing returned args
data.group.raw = pooled_f;  data.group.time = pooled_t;
data.ind.raw   = ind_f;     data.ind.time   = ind_t;
data.ind.group_id = cell_number;

infer.group.spikes = spk_pooled;   infer.group.drift = drift_pooled;    infer.group.params = parest_pooled;
infer.ind.spikes   = spk_ind;      infer.ind.drift   = drift_ind;       infer.ind.params   = parest_ind;

infer.corr.lags = lags;     infer.corr.crosscorr = lag_crosscorr;   infer.corr.smooth_bin = smooth_bin;
end

function [ lags, lag_crosscorr ] = calc_SCC( spk_pooled, spk_ind, pooled_t, ind_t, groupID, smooth_bin, dt, tmax )
% to calculate pairwise SCC between estimates of points on same
% region-of-interest
  
  smooth_factor = 2*(floor(smooth_bin*0.5/dt))+1;
  nGroups = numel(spk_pooled);
  for gg = 1:nGroups
      indx = ( groupID == gg );
%       x  = [spk_pooled{gg}{1}, cell2mat(spk_ind(indx))];  %concatenate group-and-ind traces
      x  = [spk_pooled{gg}, cell2mat(spk_ind(indx))];  %concatenate group-and-ind traces
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