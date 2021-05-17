function plot_benchmarked_traces( data, infer, varargin )

    pooled_f   = data.group.raw;        ind_f   = data.ind.raw;
    pooled_t   = data.group.time;       ind_t   = data.ind.time;
    spk_pooled = infer.group.spikes;    spk_ind = infer.ind.spikes;
    
    if nargin>3
        g = varargin{1}; r = varargin{2}; % group struct and raw data
        new_id  = []; cell_number = []; nsoma = length(g.soma); nog = length(g.other_groups);
        for ss = 1:nsoma, nid = length(g.soma{ss}(:));          new_id = [new_id;  g.soma{ss}(:) ];           cell_number = [cell_number,         ss*ones(1,nid)]; end
        for gg = 1:nog
            group = g.other_groups{gg}; 
            if iscell(group), group = group{1}; end 
            nid = length(group(:));  new_id = [new_id;  group(:) ];   cell_number = [cell_number, (nsoma+gg)*ones(1,nid)]; 
        end
        r = r(:,:, new_id); 
        nroi = length(new_id); tm = size(r,1); trl = size(r,2);

        % Find good individual ROI with high maximal dff values belonging to the
        % cells above

        zscored = (r - repmat( mean(mean(r,1),2), [tm, trl, 1] ) )./(repmat(reshape( std( reshape(r, [tm*trl, nroi]) ) ,[1 1 nroi] ),[tm,trl,1]));
        pid = prctile( reshape(zscored, [tm*trl, nroi]), 95, 1 );
        include_id = find( pid >= 0.5*prctile(pid,80) );       % Atleast half as bright as the top 20% of ROIs
        cell_number = cell_number(include_id);
    else
        cell_number = varargin{1};
    end
    
    nGroups = length(pooled_f);
    nROI    = numel(ind_f);
    toplot_traces = true;
    if nargin==5, plot_group = varargin{3}; else, plot_group = false(nGroups,1); end
    
    
    c=nROI+nGroups*3+3;
    
    if toplot_traces
        figure;hold on;
       for gg = 1:nGroups
           group_roi = find(cell_number==gg); ctr=[];
           for jj=1:length(group_roi)
                   indx = group_roi(jj); spike_raster = spk_ind{indx}; ctr = [ctr,sum(spike_raster)];
           end
           if length(group_roi)>1 && (min(ctr) > 10 || plot_group(gg) )
               gg
               c=c-3;
               % Plot spike train from averaged trace
               
               spike_raster = spk_pooled{gg}; 
               if iscell(spike_raster), spike_raster=spike_raster{1}; end
                spike_raster(spike_raster == 0) = NaN;
    %            spike_raster = spk_pooled{gg}; spike_raster(spike_raster == 0) = NaN;
               plot( pooled_t(:,gg), spike_raster+c, 'LineWidth',2 );
               c=c-1;
               % Plot spike trains from individual imaged poimt ROI
               
               for jj=1:length(group_roi)
                   indx = group_roi(jj); spike_raster = spk_ind{indx}; spike_raster(spike_raster == 0) = NaN;
        %            size(spike_raster);
                   plot( ind_t(:,indx), 0.8*double(spike_raster)+c, 'LineWidth',1 );
                   c=c-1;
               end
           end
       end
    end
end