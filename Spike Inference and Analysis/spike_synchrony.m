function [ res ] = spike_synchrony( spikes, varargin )
%%
% Adapted from SPIKY Implementation.
%
% Pre-requisites from SPIKY package:
%
%
% INPUT:
%   - spikes        [Cell Array]    Each cell is an array of spiketimes
%                                   corresponding to a single neuron.
%
%



%% Parse arguments
    p = parse_inputs( varargin{:} );
    p.num_trains = numel(spikes);      % Number of spike trains given
    p.num_spikes = cellfun( 'length', spikes );

    if ~isfield(p,'tmax'), p.tmax = max( cellfun(@(v) max([-Inf v]), spikes )); end


    % Create params structure for SPIKY-loop
    [ para.tmin, para.tmax, para.dts ] = deal( p.tmin, p.tmax, p.dt );
    m_para.all_measures_string={'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_forward';'SPIKE_synchro';'SPIKE_order';'PSTH'};  % order of select_measures
    para.select_measures  = [   p.isi   p.spikedist   p.spikedist_real   p.spikedist_forw   p.spike_synchro   p.spike_order   p.psth];
    d_para=para;    SPIKY_check_spikes
    if ret==1,      return;             end
    para=d_para;
   
    m_para.isi=1;
    m_para.spike=[2 3 4];
    m_para.spike_sync=[5 6];
    m_para.psth=7;

% plotting=7;           % +1:spikes,+2:dissimilarity profile,+4:dissimilarity matrix
% plot_profiles=3;      % 1:all,2:Groups only,3:all+groups
% 
    if any(para.select_measures)
        res = SPIKY_loop_f_distances(spikes,para)  %#ok<NOPRT>
    end


end



%% Options Parser
function params = parse_inputs( varargin )

    params.tmin             = 0;
    params.dt               = 5;%ms     % Assuming maximal acquisition rates of 200Hz
    
    
    % What measures to calculate?
    params.isi              = 0;
    params.spikedist        = 1;        % Regular SPIKE-distance calculated by default (not causal)
    params.spikedist_real   = 0;
    params.spikedist_forw   = 0;
    params.spike_synchro    = 1;        % Spike synchrony calculated by default
    params.spike_order      = 0;
    params.psth             = 0;
    
    if nargin>0
    ctr = 1;
    while ctr <= nargin
       if strcmpi( varargin{ctr}, 'tmin' )
           params.tmin = varargin{ctr+1};
           ctr = ctr+2;
       elseif strcmpi( varargin{ctr}, 'tmax' )
           params.tmax = varargin{ctr+1};
           ctr = ctr+2;
       elseif strcmpi( varargin{ctr}, 'dt' )
           params.dt = varargin{ctr+1};
           ctr = ctr+2;
       elseif strcmpi( varargin{ctr}, 'isi' )
           params.isi = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'spike' )
           params.spikedist = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'real' )
           params.spikedist_real = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'forward' )
           params.spikedist_forw = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'sync' )
           params.spike_synchro = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'order' )
           params.spike_order = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'psth' )
           params.psth = 1;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'nosync' )
           params.spike_synchro = 0;
           ctr = ctr+1;
       elseif strcmpi( varargin{ctr}, 'nospike' )
           params.spikedist = 0;
           ctr = ctr+1;
       end
       
       
    end
    end

end