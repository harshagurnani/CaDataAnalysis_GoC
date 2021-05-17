function gaptimes = get_time_bw_trials(exp_path, varargin)
%% Load speed data and locate inter-trial periods by finding where time counting restarts
% Optional trial is the upper bound for trial durations in seconds to
% remove initial junk data. Can add a different check in the future.

    %sd = dir([ exp_path, '/Speed_Data/Speed data*.txt']);
    %speed_data = importdata([ exp_path, '/Speed_Data/', sd.name]);
    sd = dir([ exp_path, 'Speed data*.txt']);
    speed_data = importdata([ exp_path, sd.name]);
    speed_data = speed_data.data(:,1:2);
    
    last = size(speed_data,1);
    flips = find(speed_data(2:last,1) - speed_data(1:last-1,1) < 0 );
    if nargin > 1
        ubound = varargin{1};
    if speed_data(flips(1)) > ubound*1000*1000
        warning('First gap greater than %d s - clipping it', ubound)
        flips = flips(2:end);
    end
    end
    n_gaps = floor(size(flips,1)/2);    %A pair of flips at each trial gap
    gaptimes = zeros(1,n_gaps);   
    
    gaptimes(1,:)= speed_data( flips(2*[1:n_gaps]))/1000;   %ms
    if max(gaptimes) > 500
       warning('Inter-trial time greater than 0.5 s')
    end
end