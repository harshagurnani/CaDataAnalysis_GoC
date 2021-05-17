function [LP] = get_loco_period( speed, abs_motion, min_LP_period, period_gap, varargin)
%% SYNTAX: Loco_period = get_loco_period( speed, abs_motion, period_gap ) OR
%                      = get_loco_period( speed, abs_motion, period_gap, mean_thresh ) OR
%                      = get_loco_period( speed, abs_motion, period_gap, mean_thresh, mean_window )
%
% Arguments:
%   - speed:            2-D Array of Time stamps and speed (raw data in Speed Data file)
%                       Can also be n_trials x 1 cell-array but will be reshaped
%                       into single array for etecting locomotion periods (because
%                       of small gap ~ 10 ms between trials)
%   - abs_motion:       (Logical) Convert speed to absolute values to calculate
%                       motion periods
%   - min_LP_period:    Minimum LP period (in ms)
%   - period_gap:       Minimum interval between 2 distinct LP periods (ms)
%
%
% {Optional args:}
%   - Varargin{1} :- mean_thresh:   2-element array - 
%                                       - Threshold of mean speed for
%                                       detecting locomotion
%                                       - Threshold for end of locomotion
%                                       period
%                                   If either value is NaN, default value
%                                   will be used for that parameter
%   - Varargin{2} :- mean_window:   Window size (in ms) for averaging speed
%                                   data  Averaged data compared against threshold 
%                                   to detect end of locomotion period
%
%
%  RETURNS:
%      LP:  3-column array with start and end times (in ms) of locomotion
%           period with mean speed in each period in the third column
%
%

    if iscell(speed)
        n_trials = size(speed,2);
        temp=[];
        
        for jj=1:n_trials
           temp = vertcat(temp,speed{jj}); 
        end
        speed = temp;  
        clear temp
    end
    
    %Default values
    def_mean_thresh = [2.8, 1.4];  %rps
    def_mean_window = 30;   %ms
    switch nargin
        case 4
            mean_thresh = def_mean_thresh;
            mean_window = def_mean_window;
        case 5
            mean_thresh = varargin{1};
            mean_window = def_mean_window;
        case 6
            mean_thresh = varargin{1};
            mean_window = varargin{2};
            
    end
    if isnan(mean_thresh(1))
        mean_thresh(1) = def_mean_thresh(1);
    end
    if isnan(mean_thresh(2))
        mean_thresh(2) = def_mean_thresh(2);
    end
    
    
    LP=[];
    
    start = 1;
    last = size(speed,1);
    if abs_motion
        speed(:,2) = abs(speed(:,2));
    end
    
    while start <= last-1
       new_LP_start = find( speed(start:last,2)>= mean_thresh(1), 1 );    %Next speed spike 
       LP_end = 0;
       curr_TS = new_LP_start + start-1;    %Actual spike index
       if isempty(curr_TS)
           start = last;
           break
       end
       
       while LP_end < 1 
           curr_time = speed(curr_TS,1);
           last_ts = find( speed(curr_TS:last,1) - curr_time > mean_window , 1) ;  %End of averaging window at current index
           if isempty(last_ts)
               LP_end = 1;  %End of current Locomotion period
               LP_start=new_LP_start+ start-1;
               LP_last = last;
               LP = vertcat( LP , [speed(LP_start,1), ...
                                   speed(LP_last ,1), ...
                                   mean(speed(LP_start:LP_last,2))] );
               start = last+1;
           else
               %mean()
               if mean( speed(curr_TS:curr_TS+last_ts-1,2)) < mean_thresh(2) 
                   LP_end = 1;  %End of current Locomotion period
                   LP_start=new_LP_start+ start-1;
                   LP_last = curr_TS+last_ts-1;
                   if speed(LP_last,1)-speed(LP_start,1) >= min_LP_period
                    LP = vertcat( LP , [speed(LP_start,1), ...
                                       speed(LP_last ,1), ...
                                       mean(speed(LP_start:LP_last,2))] );
                   end
                   start = LP_last+1;
               else
                   curr_TS = curr_TS + 1;
               end
           end
       end
    end
    
    %period_gap = 1000;  %ms
    n_periods = size(LP,1);
    if n_periods>=1
        lp_new = LP;
        p_on = LP(1,1);
        p_off = LP(1,2);
        p_num =1;
        for jj=2:size(LP,1)
            if LP(jj,1) - p_off < period_gap 
                p_off = LP(jj,2);
                if jj==n_periods
                   ind1 = find(speed(:,1)==p_on);
                   ind2 = find(speed(:,1)==p_off);
                   lp_new(p_num,:) = [p_on, p_off, mean(speed(ind1:ind2,2)) ]; 
                   p_num = p_num+1;
                end
            else
                ind1 = find(speed(:,1)==p_on);
                ind2 = find(speed(:,1)==p_off);
                lp_new(p_num,:) = [p_on, p_off, mean(speed(ind1:ind2,2)) ];
                p_on = LP(jj,1);
                p_off = LP(jj,2);
                p_num = p_num+1;
            end

        end
        LP = lp_new(1:p_num-1,:);
    end
end