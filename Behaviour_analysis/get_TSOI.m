function [ data1, time1] = get_TSOI( alldata, alltime, ptimes, extratime)
%% Crops alldata in intervals specified by ptimes (padded by extratime)
%  ptimes = start/end times of interesting periods - for eg. locomotion, whisking - and extra padding time
%   SYNTAX: [x1, t1] = get_TSOI( dff, dfftime, LocoPeriod, extrapadding )

% Interval times are considered to be in same units as alltimes (the
% timestamps corresponding to the data to be cropped).
%
% extratime can be 
%       one number (symmetrical padding on start and end of each period) 
% OR    a 2 -element array for different padding at the start and end.
% Eg. 
%       - [x1, t1] = get_TSOI( AllX, AllT, LocoPeriod, 300 ) 
% takes 300 ms extra before and after start/end points in each LocoPeriod row.
%       - [x1, t1] = get_TSOI( AllX, AllT, LocoPeriod, [500, 0] ) 
% takes 500 ms extra before each LocoPeriod row but nothing extra (0 ms) at the end.
%
% INPUT:
%
%   - alldata       [1xT vector]        Timeseries to be cropped
%   - alltimes      [1xT vector]        Timestamps for alldata
%   - ptimes        [nx2 array]         Column 1 is start time, column 2 is end time of the corresponding period of interest
%                                       Each row is a different period of interest
%   - extratime     [1 or 1x2 vector]   Extra time before and after each
%                                       POI to include in cropped series
%
% OUTPUT:
%   - data1         [nx1 cell array]    Each cell is a vector with the data
%                                       (alldata) cropped within the
%                                       correspondong period of interest
%   - time1         [nx1 cell array]    Timestamps for the correspodong
%                                       cropped data series
%-------------------------------------------------------------------------


n_periods = size(ptimes,1);
%% If only 1 element for extratime, use same for end padding.
if max(size(extratime)) == 1
    extratime= [extratime extratime];
end

% Initialise return arguments as empty
data1={};
time1={};

for prd = 1:n_periods
   start =  find(alltime > ptimes(prd,1) - extratime(1), 1);
   last = find(alltime > ptimes(prd,2) + extratime(2), 1);
   if isempty(last)
       last=max(size(alltime));
   end
   data1{prd} = alldata(start:last);
   time1{prd} = alltime(start:last);
end


end