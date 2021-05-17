function [ data, infer ] = benchmark_patch( p, r, t, masks, roi_indx, g, soma_indx, varargin )
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
nPatch = size(r, 5);
if isempty( roi_indx)
    if numel(masks) == nPatch
        roi_indx = 1:nPatch;
        disp('ROI index of masks not given: Using masks on sequential patches')
    else
        error('Cannot proceed: ROI index of masks not given, and number of masks is not equal to number of patches')
    end
end

if isempty(g),          g = num2cell(1:numel(masks)); end         % All masks are for different soma
if isempty(soma_indx),  soma_indx = 1:length(g);      end         % All masks are for soma

opts.do_infer = true;
opts.group_method = 'batch';
opts.ind_method = 'batch';
opts.do_plot = true;    opts.do_corr = true;
data_given = false;
opts.plot_what = [];
opts.corr_plot = 1;
opts.corr_dt = 100; %ms

opts_given = false;
if nargin>7
    if isstruct( varargin{1} )
        tmp = varargin{1};
        if isfield( tmp, 'group_method'), opts.group_method = tmp.group_method; end
        if isfield( tmp, 'ind_method'),   opts.ind_method   = tmp.ind_method;   end
        if isfield( tmp, 'do_infer'),     opts.do_infer     = tmp.do_infer;     end
        if isfield( tmp, 'do_plot'),      opts.do_plot      = tmp.do_plot;      end
        if isfield( tmp, 'do_corr'),      opts.do_corr      = tmp.do_corr;      end
        if isfield( tmp, 'plot_soma'),    opts.plot_what    = tmp.plot_soma;    end
        if isfield( tmp, 'corr_plot'),    opts.corr_plot    = tmp.corr_plot;    end
        if isfield( tmp, 'corr_dt'),      opts.corr_dt      = tmp.corr_dt;      end
        
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

if nargin>8
    if opts_given, data_id = 2; else, data_id = 1; end
    data = varargin{data_id}; infer = varargin{data_id+1};
    data_given = true;
end

if ~data_given && ~opts.do_infer
    opts.do_plot = false; opts.do_corr= false;
end

%% Pooling data in groups

nroi = length(g);   nsoma = length(soma_indx);
% What size sub-patches to create
soma_sizes = zeros(nsoma,1);
for soma = 1:nsoma
    for ii=g{soma_indx(soma)}
        soma_sizes(soma) = soma_sizes(soma) + sum( sum(masks{ii}) );
    end
end
subpatch = [20 40 100 150 200 250 300 400 500 750 1000];
patchrep = [ 4  4   4   3   3   3   3   3   3   3   3 ];
ii = find(subpatch>max(soma_sizes),1 ); 
if isempty(ii)||isnan(ii), ii = length(subpatch)+1; end
subpatch = subpatch(1:ii-1); patchrep = patchrep(1:ii-1);


if data_given
    subpatch = [20 40 100 150 200 250 300 400 500 750 1000];
    patchrep = [ 4  4   4   3   3   3   3   3   3   3   3 ];
    
    pooled_f = data.group.raw;
    pooled_t = data.group.time;
    
    for ii = 1:numel(data.ind)
        ind_f{ii}       = data.ind{ii}.raw;
        ind_t{ii}       = data.ind{ii}.time;
        cell_number{ii} = data.ind{ii}.group_id;
    end
    subpatch = subpatch(1:numel(data.ind));
    patchrep = patchrep(1:numel(data.ind));
    
else
    
    tm = size(r, 3); trl = size(r, 4);  %No. of timepoints and trials
    PpL = size(r,1); lns = size(r,2);   %Pixels per line, and number of lines

    % Creating full patch traces for all soma
    cropping = 0;
    [ all_traces_raw, time_traces_raw] = apply_part_mask_and_avg( r, t, cropping, masks, ...
                                            'do_avg', true, 'avg_groups', g, 'roi_indx', roi_indx, ...
                                            'do_norm', false, 'do_plot', false);
    
    pooled_f = mat2cell( reshape(all_traces_raw, tm*trl, nroi),   tm*trl, ones(1,nroi) );
    pooled_t = mat2cell( reshape(time_traces_raw, tm*trl, nroi),  tm*trl, ones(1,nroi) );

    pooled_f = pooled_f(soma_indx);
    pooled_t = pooled_t(soma_indx);
    

    % Create masks for new subgroups
    ind_f = cell(1,length(subpatch)); ind_t = ind_f;
    cell_number = ind_f; 
    pxl_number  = ind_f;

    for soma = 1:nsoma
       g_indx = soma_indx(soma);            %between 1 and length(g)
       mask_indx = g{ g_indx} ;             %mask_indx could be a vector in case of multiple patches (with its own mask) on the same soma
       patch_indx = roi_indx(mask_indx);    %which patch in raw data - could be vector
       nmasks = length(mask_indx);
       raw_tmp  = reshape(  permute( r(:,:,:,:,patch_indx),[1,5,2,3,4] ),  [nmasks*PpL, lns, tm, trl, 1]  );
       time_tmp = reshape(  permute( t(:,:,:,:,patch_indx),[1,5,2,3,4] ),  [nmasks*PpL, lns, tm, trl, 1]  );
       mask_tmp = cat( 1, masks{mask_indx} );
       for ii=1:length(subpatch)
          patch_size = subpatch(ii);    ndraws = patchrep(ii);
          if patch_size < sum(mask_tmp(:))
              small_masks = cell(1, ndraws); use_pxls = find(mask_tmp(:));
              for jj=1:ndraws
                  small_masks{jj} = false(size(mask_tmp));
                  small_masks{jj}( datasample( use_pxls, patch_size, 'Replace', false) ) = true;
              end
              [ sub_raw, sub_time] = apply_part_mask_and_avg( raw_tmp, time_tmp, cropping, small_masks, ...
                                            'do_avg', false, 'roi_indx', ones(1,ndraws), ...
                                            'do_norm', false, 'do_plot', false);
              ind_f{ii} = [ind_f{ii}, reshape(sub_raw,  tm*trl, ndraws)];
              ind_t{ii} = [ind_t{ii}, reshape(sub_time, tm*trl, ndraws)];
              cell_number{ii} = [cell_number{ii}; soma*ones(ndraws,1)];
              pxl_number{ii}  = [pxl_number{ii}; patch_size*ones(ndraws,1)];
          end
       end

    end  

    for ii=1:length(subpatch)
       ind_f{ii} = mat2cell( ind_f{ii}, tm*trl, ones(1,size(ind_f{ii},2)) ); 
    end
end


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
    run_method = opts.group_method;   % different observation noise? because different number of ROI may have been averaged.

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
    run_method = opts.ind_method;   % same observation noise --> same number of pixels of a patch are averaged

    for ii=1:length(subpatch)
        sprintf('Running spike inference on %d subpatches with %d pixels each... ZZZ ...', numel(ind_f{ii}), subpatch(ii) )
        [spk_ind{ii}, fit_ind, drift_ind{ii}, parest_ind{ii}]= deneux_spike_inference( ind_f{ii},   1/p.acquisition_rate, 'est_method', estimate, 'driftrate', driftrate_perc, 'run_method', run_method );
        % [spk fit drift parest]= deneux_spike_inference( x,
        % 1/p.acquisition_rate,  'driftrate', driftrate_perc, 'run_method', run_method, 'rise_time', risetime );
    end
    
    
end


if data_given &&  ~opts.do_infer
    [spk_pooled, drift_pooled, parest_pooled ] = deal( infer.group.spikes, infer.group.drift, infer.group.params );
    for ii=1:length(subpatch)
       [spk_ind{ii}, drift_ind{ii}, parest_ind{ii}] = deal( infer.ind{ii}.spikes, infer.ind{ii}.drift, infer.ind{ii}.params );  
    end
end

%% Plotting spike trains
if opts.do_plot
    if isempty( opts.plot_what )
        opts.plot_what = find( soma_sizes == max(soma_sizes) )';
    end
    
    nGroups = length(opts.plot_what);
    nROI    = sum(patchrep);    
    
    for gg = opts.plot_what
       c=nROI+nGroups*3+3;

       figure;hold on;
       c=c-3;
       % Plot spike train from averaged trace
       spike_raster = spk_pooled{gg}; if iscell(spike_raster), spike_raster = spike_raster{1}; end
       spike_raster(spike_raster == 0) = NaN;
       plot( pooled_t{gg}, spike_raster+c, 'LineWidth',2 );
       c=c-1;
       % Plot spike trains from all subpatches in descending order
       for ii = length(subpatch):-1:1
           group_roi = find(cell_number{ii}==gg);
           for jj=1:length(group_roi)
               indx = group_roi(jj); spike_raster = spk_ind{ii}{indx}; 
               spike_raster(spike_raster == 0) = NaN;
               plot( ind_t{ii}(:,indx), spike_raster+c, 'LineWidth',1 );
               c=c-1;
           end
           c=c-1;
       end

    end
end


%% Calculating total cross-correlation within groups
if opts.do_corr
    tmax = max(ind_t{1}(:));
    smooth_bin = opts.corr_dt;    %ms
    dt = 1000/p.acquisition_rate;  %ms
    for ii=1:length(subpatch)
        [ lags{ii}, lag_crosscorr{ii} ] = calc_SCC( spk_pooled, spk_ind{ii}, pooled_t, ind_t{ii}, cell_number{ii}, smooth_bin, dt, tmax );

    end
end


%% Plotting Spike-count correlation with patch size
if opts.do_plot
    lagzero = floor((numel(infer.corr{1}.lags)+1)/2);
    meanSCC=[];
    figure;
    
    switch opts.corr_plot
    % Style 1 - Small scatter for all pairs, and black bold for mean SCC
    % versus actual sub-patch size
      case 1
        
        for ii=1:length(subpatch)
            x=infer.corr{ii}.crosscorr(:,:,lagzero);
            x=x(x>eps & x<0.9999);

            scatter( subpatch(ii)*ones(1,numel(x)), x, 2,'MarkerFaceColor','k'); hold on
            meanSCC=[meanSCC, nanmean(x)];
        end
        plot(subpatch, meanSCC, 'ko-','LineWidth',2,'MarkerFaceColor','k')
    
    %Style 2 - SCC versus fractional subpatch size
      case 2
          all_tmpscc = []; frac_size = [];
          cmap=colormap(jet);
          for ii = 1:length(subpatch)
          for gg = 1:numel(spk_pooled)
              rows = find(cell_number{ii}==gg);
              if ~isempty(rows)
                  [ lags2, tmp_crosscorr ] = calc_SCC( spk_pooled(gg), spk_ind{ii}(rows), pooled_t(gg), ind_t{ii}(:,rows), ones(numel(rows),1), smooth_bin, dt, tmax );
                  lagzero = floor((numel(lags2)+1)/2);
                  tmp_crosscorr = tmp_crosscorr(:,:,lagzero);
                  tmp_crosscorr = tmp_crosscorr(tmp_crosscorr > 0 & tmp_crosscorr < 0.9999 );
                  all_tmpscc = [all_tmpscc; tmp_crosscorr ];
                  frac_size = [frac_size; ones(numel(tmp_crosscorr),1)* subpatch(ii)/soma_sizes(gg) ];
                  scatter( ones(numel(tmp_crosscorr),1)* subpatch(ii)/soma_sizes(gg), tmp_crosscorr, 2, 'MarkerFaceColor', cmap(gg,:) ); hold on
              end
          end
          end
          
          bins = 0:0.1:(max(frac_size)+0.1);
          for jj = 2:numel(bins)
              x = all_tmpscc( frac_size>=bins(jj-1) & frac_size <bins(jj) );
              meanSCC = [meanSCC, nanmean(x) ];
          end
          bins = bins(2:end);
          plot( bins(~isnan(meanSCC))-0.05, meanSCC(~isnan(meanSCC)), 'ko-' ) 
        
            
    end    
end



%% Preparing returned args

data.group.raw = pooled_f;  data.group.time = pooled_t;
for ii=1:length(subpatch)
    data.ind{ii}.raw   = ind_f{ii};             data.ind{ii}.time   = ind_t{ii};
    data.ind{ii}.group_id = cell_number{ii};    %data.ind{ii}.patch_size = pxl_number{ii};
end

if opts.do_infer
    infer.group.spikes = spk_pooled;   infer.group.drift = drift_pooled;    infer.group.params = parest_pooled;
    for ii = 1:length(subpatch)
        infer.ind{ii}.spikes   = spk_ind{ii};      infer.ind{ii}.drift   = drift_ind{ii};       infer.ind{ii}.params   = parest_ind{ii};
    end
end

if opts.do_corr
    for ii=1:length(subpatch)
        infer.corr{ii}.lags = lags{ii};     infer.corr{ii}.crosscorr = lag_crosscorr{ii};   infer.corr{ii}.smooth_bin{ii} = smooth_bin;
    end
end
% 
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
      xt = [pooled_t{gg}, ind_t(:,indx)];
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