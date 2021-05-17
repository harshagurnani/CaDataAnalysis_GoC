function [ speed, TrialTime ] = read_speed_labview( fullpath, varargin )
%% speed = get_speed_data( 'path\speed_file.txt') OR 
%        = get_speed_data( 'path\speed_file.txt', T ) if you want to check that first gap doesn't come after more than 'T' seconds. 

    if nargin>1
        lengthchck = varargin{1}; %upper limit on length of first trial (in sec) - to exclude rubbish data in the beginning.
    end
    
    
    if length(fullpath)>4 && strcmp(fullpath(length(fullpath)-3:end),'.txt')
            speed_data = importdata(fullpath);
    else
        sd = dir([ fullpath, '/Speed data*.txt']);
        speed_data = importdata([ fullpath, '/',sd.name]);
    end
    
    absTime = speed_data.textdata(:,2); 
    relTime = speed_data.data(:,1)*.001; %in ms

    % convert absTime from string to num
    absTime = arrayfun( @(jj) ( str2double(absTime{jj}(1:2))*60*60 + str2double(absTime{jj}(4:5))*60 + str2double(absTime{jj}(7:end)) )*1000, 1:length(absTime));% in ms
    absTime = absTime';

    % find encoder resets - start and end of trial
    flips = find(relTime(2:end)<relTime(1:end-1))+1;    

    exp_start = 1;  % experimental start trigger
    % check for delayed trigger?
    tlen = prctile(diff(flips), 80 ); 
    if relTime(flips(1)-1) < 0.5*tlen, exp_start = flips(1); flips = flips(2:end);   end
    absTime = absTime - absTime(exp_start) + relTime(flips(1)); %set exp_start trigger time as absolute zero 

    nTrials = floor(length(flips)/2)+1;
    % Add 2.3ms per reset
    if nargin==3,   reset_shift = varargin{3};
    else,           reset_shift = 2.3; %ms
    end
    for jj=2:nTrials
        absTime( flips(2*(jj-1)) : end ) = absTime( flips(2*(jj-1)) : end ) + reset_shift;
    end
    
    % Trial start and end time - absolute times
    TrialTime = nan( nTrials, 2 );
    TrialTime(1,1) = 0; TrialTime(1,2) = absTime( flips(1) ); 
    for jj=2:nTrials
        TrialTime(jj,1) = absTime(flips(2*(jj-1)))-relTime(flips(2*(jj-1))); % Trial start time = abstime at flip - rel time at flip
        if jj==nTrials, TrialTime(jj,2) = absTime(end); 
        else, TrialTime(jj,2) = absTime(flips(2*jj-1));                      % trial end time = abstime at next flip
        end
    end

    speed_data = speed_data.data(:,1:2);
    
    last = size(speed_data,1);
    
    if nargin==2
        if speed_data(flips(1)) > lengthchck*1000*1000
        warning('First gap greater than %s s - clipping it', lengthchck)
        
        speed_data=speed_data(flips(1)+1:end,1:2);
        flips = flips(2:end)-flips(1);
        end
        last = size(speed_data,1);
    end
    n_flips = size(flips,1);
    n_gaps = floor(n_flips/2);    %A pair of flips at each trial gap
    
    %%%No dead period after last trial?
    if mod(n_flips,2) == 0
        flips = [flips;last];
    end
    gaptimes = zeros(1,n_gaps);   
    
    prev_trial_length = speed_data(flips(2*[0:n_gaps]+1),1); %us
    gaptimes(1,:)= speed_data( flips(2*[1:n_gaps]))-2000;   %us
    
    speed = {};
    speed{1} = speed_data(1:flips(1),:);    %us
    speed{1}(:,1)=speed{1}(:,1)/1000;       %ms
    for jj=1:n_gaps
       total_prev_time = sum(prev_trial_length(1:jj)) ... Previous trial length
                            + sum( gaptimes(1:jj))  ;   % Gap between trials
       speed{jj+1} = speed_data(flips(2*jj)+1:flips(2*jj+1),:) + ...
                        [total_prev_time,0];    %us
       speed{jj+1}(:,1)=speed{jj+1}(:,1)/1000;  %ms
    end
end