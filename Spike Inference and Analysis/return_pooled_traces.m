function varargout = return_pooled_traces( data, infer, varargin )
%% To return or plot traces pooled across pixels from same cell  - based on cell ids of individual pixels, or recovering this list from group structure and raw traces 
% If no return argument is needed, it pools and plots spike trains for
% groups and individual pixels. Otherwise it returns 'traces' - reordered
% spike trains - the group averaged inference followed by individual
% pixels. By default, it does not plot if return argument is asked, but
% additional options can be used to plot as well.
%
%   Usage:
%   ------------
%       1) STRONGLY RECOMMENDED !!!!!
%          [(traces (, group_id ] = return_pooled_traces( data, infer );
%                        Read cell IDs of individual pixels in data.ind.group_id
%       2) [(traces (, group_id ] = return_pooled_traces( data, infer, g, raw );
%                        Use structure g to group pixels, zscore data to
%                        identify usable pixels to get roi corresponding to
%                        inferred spike trains in data/infer.ind. Criteria
%                        to select points using their maximal zscores
%                        should be same as the original criteria (to keep
%                        roi in data/infer). 
%
%   Inputs:
%   ------------
%       - data      As returned by function < benchmark_multiple_roi >
%       - infer     As returned by function < benchmark_multiple_roi >
%
%   (Optional): 
%   ------------
%   Must be given as name-value pairs
%       - 'do_norm'     Boolean     To return pooled DFF (as fourth output ONLY)
%       - 'do_plot'     Boolean     If you want to plot pooled traces.
%                                   Default: True if there are no return outputs, o/w false
%       - 'raster'      Boolean     If spike trains to be plotted. True by default.
%                                   Only used if 'do_plot' is true.
%       - 'dff'         Boolean     If plotting normalised traces instead!
%                                   (turns raster OFF)
%       - 'scale'       Scalar      Y-axis scale (spacing between traces) -
%                                   Default: 1 for spike trains, 0.2 for
%                                   dff.
%       - 'group_size'  Scalar      Minimum number of ind pixels in ROI to
%                                   keep it for plotting
%       - 'plot_nrn'    1xN Vector  Neuron/ROI IDs to be kept for plotting                            
%
%
%
%   Outputs (Optional):
%   ------------
%       1] traces       Cell array  Each cell has spike-inference results
%                                   from a single averaged/individual pixel
%                                   data. Could be the MAP estimate - a
%                                   single spike train, or multiple
%                                   samples. The traces are ordered such
%                                   that the inference of the
%                                   group-averaged data is followed by
%                                   inference from individual pixels.
%       2] group_id     Array       Somatic (ROI) ID of each cell in 'traces'.
%                                   The first occurrence of any somatic id
%                                   is the group-average data of that ROI,
%                                   and then all picels that make that ROI.
%       3] times        Matrix      Acquisition timepoints for each
%                                   group/pixel in traces
%       4] traces_dff   Matrix      Raw fluorescence OR Normalized and
%                                   smoothed DFF - of pooled data in same
%                                   order as traces and times
%
%   REFERENCE NOTES:
%   --------------------
%   The structure g has fields: 'soma' and 'other_groups' ('dendrite' field
%   is not used here)
%   Each field is a cell array - each cell consists a list (array) of
%   pixel numbers corrsponding to the ROI. 
%
%
%
%   Harsha G.
%   Nov 2017

    pooled_f   = data.group.raw;	pooled_t   = data.group.time;	spk_pooled = infer.group.spikes; % group-averaged data
    ind_f   = data.ind.raw;     	ind_t   = data.ind.time;     	spk_ind = infer.ind.spikes;      % individual pixel data
        
    optParams.do_plot = (nargout==0);
    [g, r, optParams] = parse_args( optParams, varargin{:} );
    if ~isvector(   spk_ind{1} ), optParams.do_plot = false; disp('Spike inference did not return a single spike train - will not generate plot'); end
    
    %% Re-evaluate ROI ID of inferred spike trains
    if ~isempty(g)
        % Have to find good roi from scratch and find cell_numbers of
        % included roi using group structure 'g' and raw data 'r'
        
        new_id  = []; cell_number = []; nsoma = length(g.soma); nog = length(g.other_groups);
        for ss = 1:nsoma
            nid = length(g.soma{ss}(:));    new_id = [new_id;  g.soma{ss}(:) ];     cell_number = [cell_number,      ss   *ones(1,nid)]; 
        end
        for gg = 1:nog
            group = g.other_groups{gg};     if iscell(group), group = group{1}; end 
            nid = length(group(:));         new_id = [new_id;  group(:) ];          cell_number = [cell_number, (nsoma+gg)*ones(1,nid)]; 
        end
        
        r = r(:,:, new_id);     tm = size(r,1); trl = size(r,2);
        nroi = length(new_id); 

        % Find good individual ROI with high maximal dff values belonging to the
        % cells above

        zscored = (r - repmat( mean(mean(r,1),2), [tm, trl, 1] ) )./(repmat(reshape( std( reshape(r, [tm*trl, nroi]) ) ,[1 1 nroi] ),[tm,trl,1]));
        pid = prctile( reshape(zscored, [tm*trl, nroi]), 95, 1 );
        include_id = ( pid >= 0.5*prctile(pid,80) );       % Atleast half as bright as the top 20% of ROIs
        cell_number = cell_number(include_id);
    else
        cell_number = data.ind.group_id;
    end
    
    %% Pool (reorder) Traces
    
    nGroups  = length(pooled_f);
    nPixROI  = numel(ind_f);
     
    
    traces = cell(1, nGroups+nPixROI);
    group_id = nan( 1, nGroups+nPixROI );
    times = nan( size(pooled_t,1), nGroups+nPixROI ); 
    if (~optParams.raster && optParams.do_plot) || (nargout>3),      poolraw = true; traces_dff = traces;        
    else,   poolraw = false;    end
    if (~optParams.raster && optParams.do_plot),	optParams.do_norm = true; optParams.dff = true; end
    
    ctr = 1; 
    for gg = 1:nGroups
       
       % spike train from group-averaged data
       spike_raster = spk_pooled{gg}; 
       if iscell(spike_raster), spike_raster=spike_raster{1}; end
       traces{ctr}  = spike_raster;
       times(:,ctr) = pooled_t(:, gg);
       if poolraw,  traces_dff{ctr} = data.group.raw{gg}; end
       group_id(ctr) = gg;  ctr = ctr+1;
       
       % spike trains from individual pixels
       group_roi = find(cell_number==gg); 
       for jj=1:length(group_roi)
           indx = group_roi(jj); spike_raster = spk_ind{indx}; 
           traces{ctr}  = spike_raster; 
           times(:,ctr) = ind_t(:, indx);
           if poolraw,  traces_dff{ctr} = data.ind.raw{indx}; end
           group_id(ctr) = gg;  ctr = ctr+1; 
       end
    end
    if poolraw,             traces_dff      = cell2mat(traces_dff);  end
    if optParams.do_norm,   tm = size(traces_dff,1); rn = size(traces_dff,2);
                            [traces_dff, ~] = normalize_and_smooth( reshape(traces_dff,[tm,1,rn]), reshape(times,[tm,1,rn]), [] );
                            traces_dff = squeeze(traces_dff);
    end
    
    if nargout >0,	varargout{1} = traces;      end
    if nargout >1,  varargout{2} = group_id;    end
    if nargout >2,  varargout{3} = times;       end
    if nargout >3,  varargout{4} = traces_dff;  end
    
    
    %% Plot Traces?
    if optParams.do_plot
        figure();hold on; 
        c =  (nPixROI+nGroups*3+3) ;

        if isfield( optParams, 'group_size' )
            % only plot groups with minimum 'group_size' number of individual
            % pixels
            optParams.plot_nrn = [];
            for gg=1:nGroups
                if sum(group_id==gg) > optParams.group_size, optParams.plot_nrn = [optParams.plot_nrn, gg]; end
            end
        elseif ~isfield( optParams, 'plot_nrn' )
            % if no nrn/group ids are given, plot all groups
            optParams.plot_nrn = 1:nGroups;
        end

        for gg = optParams.plot_nrn
            all_id = find( group_id == gg );
            group_avg = all_id(1);
            ind_pixel = all_id(2:end); 

            % Plot spike train from group-averaged fluorescence
            if optParams.raster,   plot_trace = traces{group_avg} ;
            elseif optParams.dff,  plot_trace = traces_dff(:,group_avg);     end
            plot( times(:,gg), plot_trace+ optParams.scale*c,                   'LineWidth',2 );
            c=c-1; 

            % Plot spike trains from individual imaged point ROI        
            for jj=ind_pixel
                if optParams.raster,   plot_trace = traces{jj} ;
                elseif optParams.dff,  plot_trace = traces_dff(:,jj);     end
                plot( times(:,jj), 0.8*double(plot_trace)+ optParams.scale*c,   'LineWidth',1 );
                c=c-1;
            end

            c = c-2;
        end
    end

end

function [ g, raw, optParams] = parse_args( optParams, varargin )

    g=[]; raw = [];
    optParams.do_norm = false; % Normalize pooled Fluorescence data?
    optParams.raster = true;   % A single spike train in inference results - spikes in each time bin
    optParams.dff = false;     % Plot DFF of group-averages and ind pixels?
    
    if ~isempty(varargin)
       if isstruct( varargin{1} ),  g = varargin{1}; raw = varargin{2}; read_var =  3;
       else,                        read_var = 1;       end
       
       while read_var <= numel(varargin)
            optParams.(varargin{read_var}) = varargin{read_var+1};
            read_var = read_var + 2;
       end
    end
    
    % make dff as toggle
    if optParams.dff, optParams.raster = false;  end
    
  
    if ~isfield( optParams, 'scale' )
       if optParams.raster, optParams.scale = 1;
       else,                optParams.scale = 0.4;
       end
    end
end