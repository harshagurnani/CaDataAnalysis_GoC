function [ speed ] = get_speed_rmvsample( fullpath, sample, varargin )
%% speed = get_speed_data( 'path\speed_file.txt') OR 
%        = get_speed_data( 'path\speed_file.txt', T ) if you want to check that first gap doesn't come after more than 'T' seconds. 

    if ~isempty(varargin)
        lengthchck = varargin{1}; %upper limit on length of first trial (in sec) - to exclude rubbish data in the beginning.
    end
    
    if length(fullpath)>4 && strcmp(fullpath(length(fullpath)-3:end),'.txt')
            speed_data = importdata(fullpath);
    else
        sd = dir([ fullpath, '/Speed data*.txt']);
        speed_data = importdata([ fullpath, '/',sd.name]);
    end
    speed_data = speed_data.data(:,[1 3]);
    
    last = size(speed_data,1);
    flips = find(speed_data(2:last,1) - speed_data(1:last-1,1) < 0 );
    
    if ~isempty(varargin)
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
    gaptimes(1,:)= speed_data( flips(2*[1:n_gaps]));   %us
    
    %undersampling factor
    speed = cell(1,n_gaps+1);
    speed{1} = speed_data(1:sample:flips(1),:);    %us
    speed{1}(:,1)=speed{1}(:,1)/1000;       %ms
    for jj=1:n_gaps
       total_prev_time = sum(prev_trial_length(1:jj)) ... Previous trial length
                            + sum( gaptimes(1:jj))  ;   % Gap between trials
       speed{jj+1} = speed_data(flips(2*jj)+1:sample:flips(2*jj+1),:) + ...
                        [total_prev_time,0];    %us
       speed{jj+1}(:,1)=speed{jj+1}(:,1)/1000;  %ms
    end
    
    % At this point, it has angles for all trials (and has been downsampled
    % by sample). Now compute differences
    for jj=1:n_gaps+1
       tmp =  diff(speed{jj}(:,2));
       x1 = ( abs(tmp) < abs(tmp+360) );
       x2 = (abs(tmp) < abs(tmp-360) );
       tmp = (tmp.*(x1&x2) + (tmp+360).*(~x1) + (tmp-360).*(~x2));
       speed{jj}(2:end,2) = tmp;
       if jj==1, speed{jj}(1,2)= 0; else, speed{jj}(1,2)=speed{jj-1}(end,2); end
    end
end