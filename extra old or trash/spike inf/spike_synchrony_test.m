function [ res ] = spike_synchrony_test( spikes, varargin )
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

    if ~isfield(p, 'tmax'), p.tmax = max( cellfun(@(v) max([-Inf v]), spikes )); end

%% Prepare spike trains
    warnon = true;
    for jj=1:p.num_trains    
        if length(unique(spikes{jj}))  <  p.num_spikes(jj)
            if warnon, warning ('Multiple spikes with same time found - Will be collapsed into a single spike time '); warnon = false; end
        end
%         while length(unique(spikes{jj}))  <  p.num_spikes(jj)
%             % Multiple spikes with same time - jitter duplicated times by a
%             % millisecond or less 
%             if warnon, warning ('Multiple spikes with same time found - jittering the duplicated times by 0.2-1 ms'); warnon = false; end
%             [~, utime, ~ ] =unique(spikes{jj});        dup_ind = setdiff(1:p.num_spikes, utime );
%             jitter = 2*rand(1,length(dup_ind))-1;      jitter = sign(jitter)*0.2 + 0.8*jitter;
%             spikes{jj}(dup_ind) = spikes{jj}(dup_ind) + jitter;
%         end
        spikes{jj} = sort(spikes{jj});      % Ensure it's still sorted.
    end

    %-------
    uspikes=cell(1,p.num_trains);      % Unique spikes
    for trac=1:p.num_trains
        uspikes{trac}=spikes{trac}(spikes{trac}>=p.tmin-100*eps & spikes{trac}<=p.tmax+100*eps);  %Crop between tmin and tmax
        uspikes{trac}=unique([p.tmin uspikes{trac} p.tmax]);      % Append tmin and tmax
    end
    num_uspikes=cellfun('length',uspikes);
    max_num_uspikes=max(num_uspikes);

    %-------
    % Real first and late spikes close to tmin and tmax? (checking spikes and
    % not uspikes) 
    non_empties=find(cellfun('length',spikes)>0);
    tmin_spikes=non_empties(abs(cellfun(@(v) v(1), spikes(cellfun('length',spikes)>0))-p.tmin)<1e-20);
    tmax_spikes=non_empties(abs(cellfun(@(v) v(end), spikes(cellfun('length',spikes)>0))-p.tmax)<1e-20);
    m_res.num_tmin_spikes=length(tmin_spikes);
    m_res.num_tmax_spikes=length(tmax_spikes);

    % Combine all spiketimes, sort and get the associated cell numbers for each
    % time
    dummy=[0 num_uspikes];
    all_indy=zeros(1,sum(num_uspikes));
    for trac=1:p.num_trains
        all_indy(sum(dummy(1:trac))+(1:num_uspikes(trac)))=trac*ones(1,num_uspikes(trac));
    end
    [all_spikes,indy]=sort([uspikes{:}]);
    num_all_spikes=length(all_spikes);

    all_trains=all_indy(indy);
    all_trains(setdiff(1:p.num_trains,tmin_spikes))=0;                         %Tmin/Tmax spikes that were appended should get cell ID = 0
    all_trains(end-p.num_trains+setdiff(1:p.num_trains,tmax_spikes))=0;

    % Also: Keep one Tmin/Tmax spike to the pooled population spiketimes if
    % there are no spikes close to Tmin and Tmax.
    if m_res.num_tmin_spikes==0 && m_res.num_tmax_spikes==0
        all_indy=p.num_trains:num_all_spikes-p.num_trains+1;      % Remove all but one appended spikes at tmin and tmax
    elseif m_res.num_tmin_spikes==0
        all_indy=[p.num_trains:num_all_spikes-p.num_trains num_all_spikes-p.num_trains+tmax_spikes];
    elseif m_res.num_tmax_spikes==0
        all_indy=[tmin_spikes p.num_trains+1:num_all_spikes-p.num_trains+1];
    else
        all_indy=[tmin_spikes p.num_trains+1:num_all_spikes-p.num_trains num_all_spikes-p.num_trains+tmax_spikes];
    end

    all_spikes=all_spikes(all_indy);
    all_trains=all_trains(all_indy);
    
    

%% Calculating isi
    m_res.all_isi           = diff(all_spikes);                 % All interspike intervals
    m_res.num_all_isi       = length(all_spikes)-1;
    
    m_res.isi               = m_res.all_isi(m_res.all_isi>0);   % Nonzero ISI
    m_res.num_isi           = length(m_res.isi);
    
    [ m_res.cum_isi, indx, ~] = unique(all_spikes);           %All times to compute at
    
    clear all_indy indy 

    times = m_res.cum_isi;%sort(unique([m_res.cum_isi, p.tmin:(p.tmax-p.tmin)/10000:p.tmax] ));
    res.SpikeDist_times = times;
    nT = length(times);
    
    [ TPrev, TForw, xPrev, xForw ]  = deal( nan( nT-1, p.num_trains ) );  
    [ delTPrev, delTForw ]          = deal( nan( nT-1, p.num_trains, p.num_trains ) );  % Bi-variate premeasure
    [ delTPrev2, delTForw2 ]        = deal( nan( nT-1, p.num_trains, p.num_trains ) );  % Premeasures for improved Spike-distance
      delTPrev_real                 = nan( nT-1, p.num_trains, p.num_trains ) ;         % Premeasures for realtime Spike-distance
    
    [ SDist1, SDist2, SDist, SDist_real ] = deal( nan( nT-1, p.num_trains, p.num_trains ) );  % Spike-distance measures
    
    % These variables only change value at spiketimes and need only be
    % evaluated there. Complexity is O(N) where N is number of total spikes
    for tt = 1:nT-1
        TPrev( tt, :) = cellfun( @(v) find( v<=times(tt), 1, 'last') , uspikes );      % Last previous spike before time tt
        TForw( tt, :) = cellfun( @(v) find( v> times(tt), 1)         , uspikes );      % Following spike
        
        xPrev( tt, :) = times(tt) - TPrev( tt, :);                      % Local weights
        xForw( tt, :) = TForw( tt, :) - times(tt);
        
    end
    xISI = TForw - TPrev;     % ISI of each cell

%% Instantaneous spiketime distance 
    % Normalised by instantaneous isi within each train

    if p.SpikeDist.old  || p.SpikeDist.local
        % old or locally-weighted spike distance
        for trac1 = 1:p.num_trains
            for trac2 = trac1+1 : p.num_trains
                % Absolute isi
                delTPrev( :, trac1, trac2) = abs( TPrev(:, trac1) - TPrev(:, trac2) );
                delTPrev( :, trac2, trac1) = delTPrev( :, trac1, trac2);
                delTForw( :, trac1, trac2) = abs( TForw(:, trac1) - TForw(:, trac2) );
                delTForw( :, trac2, trac1) = delTForw( :, trac1, trac2);

                % Old measure is a simple average
                SDist1( :, trac1, trac2 ) = squeeze( delTPrev(:, trac1, trac2) + delTForw(:, trac1, trac2) )./squeeze( xISI(:, trac1) + xISI(:,trac2) );
                SDist1( :, trac2, trac1 ) = SDist1( :, trac1, trac2 );

                % SDist2 is refined to be local weighted average
                SDist2( :, trac1, trac2 ) = 2*( delTPrev(:,trac1,trac2).*(xForw(:,trac1) + xForw(:,trac2))  + delTForw(:,trac1,trac2).*(xPrev(:,trac1) + xPrev(:,trac2)) )  ./ ...
                                                                    ( xISI(:, trac1) + xISI(:,trac2) ).^2 ;
                SDist2( :, trac2, trac1 ) = SDist2( :, trac1, trac2 ); 
            end
        end
        
        res.SpikeDist_old.matrix   = SDist1;
        res.SpikeDist_local.matrix = SDist2;
        
    end

    if p.SpikeDist.improved
        for trac1 = 1:p.num_trains
            for trac2 = trac1+1 : p.num_trains
                % for improved spike distance - compare nearest neighbour,
                % irrespective of whether it is previous or following
                delTPrev2( :, trac1, trac2 ) = arrayfun( @(v, v1, v2) min( abs(v-v1), abs(v-v2) ), TPrev(:, trac1), TPrev(:, trac2), TForw(:, trac2) );
                delTPrev2( :, trac2, trac1 ) = arrayfun( @(v, v1, v2) min( abs(v-v1), abs(v-v2) ), TPrev(:, trac2), TPrev(:, trac1), TForw(:, trac1) );
                delTForw2( :, trac1, trac2 ) = arrayfun( @(v, v1, v2) min( abs(v-v1), abs(v-v2) ), TForw(:, trac1), TPrev(:, trac2), TForw(:, trac2) );
                delTForw2( :, trac2, trac1 ) = arrayfun( @(v, v1, v2) min( abs(v-v1), abs(v-v2) ), TForw(:, trac2), TPrev(:, trac1), TForw(:, trac1) );    
                
                % SDist is the 'improved version' and used by default
                % see Kreuz et al 2013 for definition
                % each train is normalised separately before combining for
                % computing bivariate distance
                temp1 = (  delTPrev2( :, trac1).*xForw(:, trac1)  +  delTForw2( :, trac1).*xPrev(:, trac1)  ) ./ xISI(:, trac1);    
                temp2 = (  delTPrev2( :, trac2).*xForw(:, trac2)  +  delTForw2( :, trac2).*xPrev(:, trac2)  ) ./ xISI(:, trac2);    

                SDist(:, trac1, trac2) = 2*( temp1.*xISI(trac2) + temp2 .*xISI(trac1) )  ./  ( xISI(trac1) + xISI(trac2) ).^2;
                SDist(:, trac2, trac1) = SDist(:, trac1, trac2);
                
            end
        end
        res.SpikeDist.matrix = SDist;
    end
            
    if p.SpikeDist.realtime        
        % for realtime measure
        for trac1 = 1:p.num_trains
            for trac2 = trac1+1 : p.num_trains
                                                             
            for tt = 1:length( m_res.cum_isi )-1
               delTPrev_real( tt, trac1, trac2 ) = min( abs( TPrev(tt,trac1) - (uspikes{trac2}<m_res.cum_isi(tt)) ) ); 
               delTPrev_real( tt, trac2, trac1 ) = min( abs( TPrev(tt,trac2) - (uspikes{trac1}<m_res.cum_isi(tt)) ) ); 
            end
            
            % Realtime spike distance
            SDist_real(:, trac1, trac2) = 0.5*(  delTPrev_real(:, trac1, trac2) + delTPrev_real(:, trac2, trac1)  ) ./ ( xPrev(:,trac1) + xPrev(:,trac2) ) ;
            SDist_real(:, trac2, trac1) = SDist_real(:, trac1, trac2);
            end
        end
        
        res.SpikeDist_real.matrix = SDist_real;
    end
    
   
    
    for tt = 1:nT-1
%         spiking_cells = all_trains( indx(tt):indx(tt+1)-1 );
        [tmp1, tmp2, tmp, tmp_real ] = deal( [] );
        for trac1 = 1:p.num_trains%setdiff(1:p.num_trains, spiking_cells )%1:p.num_trains%
            for trac2 = trac1+1:p.num_trains%setdiff( 1:p.num_trains, [trac1 spiking_cells] )%1:p.num_trains%
               if p.SpikeDist.old,          tmp1        = [tmp1     res.SpikeDist_old.matrix(tt, trac1, trac2) ];       end 
               if p.SpikeDist.local,        tmp2        = [tmp2     res.SpikeDist_local.matrix(tt, trac1, trac2) ];     end 
               if p.SpikeDist.improved,     tmp         = [tmp      res.SpikeDist.matrix(tt, trac1, trac2) ];           end 
               if p.SpikeDist.realtime,     tmp_real	= [tmp_real res.SpikeDist_real.matrix(tt, trac1, trac2) ];      end 
            end
        end
        if p.SpikeDist.old,         res.SpikeDist_old.profile(tt)   = nanmean( tmp1 );      end 
        if p.SpikeDist.local,       res.SpikeDist_local.profile(tt) = nanmean( tmp2 );      end 
        if p.SpikeDist.improved,    res.SpikeDist.profile(tt)       = nanmean( tmp );       end 
        if p.SpikeDist.realtime,	res.SpikeDist_real.profile(tt)  = nanmean( tmp_real );	end 
        
    end
%% Spike-Synchrony


end

function params = parse_inputs( varargin )

    params.tmin = 0;
    params.num_runs = 1;
    
    params.SpikeDist.old        = true;
    params.SpikeDist.local      = true;
    params.SpikeDist.improved   = true;
    params.SpikeDist.realtime   = true;
    if nargin>0
    ctr = 1;
    while ctr <= nargin
       if strcmp( varargin{ctr}, 'tmin' )
           params.tmin = varargin{ctr+1};
           ctr = ctr+2;
       elseif strcmp( varargin{ctr}, 'tmax' )
           params.tmax = varargin{ctr+1};
           ctr = ctr+2;
       elseif strcmp( varargin{ctr}, 'nruns' )
           params.num_runs = ceil(max(varargin{ctr+1},1));
           ctr = ctr+2;


       end
    end
    end

end