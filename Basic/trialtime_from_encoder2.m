function [TrialTime, gaptimes ] = trialtime_from_encoder( params, varargin )
%% read times in encoder data, get gintertrial dur, and start and end times of trials

    if nargin>1
        lengthchck = varargin{1}; %upper limit on length of first trial (in sec) - to exclude rubbish data in the beginning.
    end
    sd = dir([ params.exp_path, '/Speed_Data/Speed data*.txt']);    %to allow for mild name changes
    speed_data = importdata([  params.exp_path, '/Speed_Data/', sd.name]);
    
    absTime = speed_data.textdata(:,2); 
    relTime = speed_data.data(:,1)*.001; %in ms

    % convert absTime from string to num
    absTime = arrayfun( @(jj) ( str2double(absTime{jj}(1:2))*60*60 + str2double(absTime{jj}(4:5))*60 + str2double(absTime{jj}(7:end)) )*1000, 1:length(absTime));% in ms
    absTime = absTime';

    % find encoder resets - start and end of trial
    flips = find(relTime(2:end)<relTime(1:end-1))+1;    

    exp_start = 1;  % experimental start trigger
    if ~isempty( flips )
    % check for delayed trigger?
        tlen = prctile(diff(flips), 80 ); 
        if relTime(flips(1)-1) < 0.5*tlen, exp_start = flips(1); flips = flips(2:end);   end
    else
        % no trigger found - only 1 trial(?)
        flips = exp_start;
    end
    absTime = absTime - absTime(exp_start) + relTime(flips(1)); %set exp_start trigger time as absolute zero 

    nTrials = floor(length(flips)/2)+1;
    if params.maxTrial > nTrials, warning('Not enough resets found in encoder data'); end
    
    % Add 2.3ms (or preset value) per encoder reset 
    try
        reset_shift = params.encoder_shift_dt;
    catch
        reset_shift = 2.3;  %ms
    end
    for jj=2:nTrials
        absTime( flips(2*(jj-1)) : end ) = absTime( flips(2*(jj-1)) : end ) + reset_shift;  %shift the subsequent trials
        relTime( flips(2*(jj-1))-1 ) = relTime( flips(2*(jj-1))-1 ) + reset_shift;      %add the missing time to intertrial end
    end
    
    % Trial start and end time - absolute times
    TrialTime = nan( nTrials, 2 );
    TrialTime(1,1) = 0; TrialTime(1,2) = absTime( flips(1) ); 
    for jj=2:nTrials
        TrialTime(jj,1) = absTime(flips(2*(jj-1)))-relTime(flips(2*(jj-1))); % Trial start time = abstime at flip - rel time at flip
        if (jj==nTrials && length(flips)<(2*nTrials-1)), TrialTime(jj,2) = absTime(end); 
        else, TrialTime(jj,2) = absTime(flips(2*jj-1));                      % trial end time = abstime at next flip
        end
    end
    
    n_flips = size(flips,1);
    n_gaps = floor(n_flips/2);    %A pair of flips at each trial gap
    
    gaptimes = zeros(n_gaps, 1);   
    gaptimes(:,1)= relTime( flips(2*[1:n_gaps])-1 );   %ms
    

end
