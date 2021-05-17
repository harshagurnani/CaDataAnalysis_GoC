function [A1, T1] = resample_ds( data1, TS1, DS_rate, interval1 )
%%                    
% ARGUMENTS:
%   - data1   [n x 1 (FLOAT)]:    
%           First data vector. Could be dFF  or speed or whisker position
%   - TS1     [n x 1 (FLOAT)]:       
%           Time stamps for first vector (ms)
%   - DS_rate [INT]:
%           Rate of downsampled data for both vectors. (Hz) 
%   - interval1 [1x2]:
%           Interval to be sampled in (start/stop time) (ms)
%   optional:
%
%
%   OUTPUT:
%
%       - A1 [1-D array]
%           Resampled array
%       - T1 [1-D array]
%           Resampled timepoints
%
% Harsha Gurnani. March 2017

    
    %% Interpolate data
    isi1 = ceil( mean(TS1(2:end)-TS1(1:end-1) ));       %F
    %%%%%%%%% Ceil or floor?
    
    %Downsampling factors need to be integer (in decimate function) so
    %resample at one of the divisor ISIs
    isi_to_use = divisors(1000/DS_rate);        % Allowed ISI [all divisors of 1000 ms]
    if (isi1 > 2*isi_to_use(end)  )
        error('Downsampling rate is too high')
    end
    
    %Closest isi to use for resampling the vectors before decimation
    new_isi = find( isi_to_use>=isi1,1);
    if isempty(new_isi)
        isi1 = isi_to_use(end);
    else
        isi1=isi_to_use(new_isi);
    end
    
    if isempty(interval1), interval1 = [TS1(1) TS1(end)]; end
    T1 = [interval1(1):isi1:interval1(2)];
    %Resample to one of the 'allowed ISI' with spline interpolation
    % Linear interpolation is better if you don't want to introduce negative values in
    % something that is supposed to be positive eg. spike counts which are
    % integral valued
    interp_method = 'linear';
    if any( mod(data1(2:end)-data1(1:end-1), 1 ) ~= 0 ) & ~any(data1<0) 
    %Non-integral jumps in data and negative values ---> unlikely to be spikes/firing rate ---> use spline
        interp_method = 'spline';
    end
    new_data1 = interp1( TS1, data1,  T1 , interp_method);     
    
    
    %% Filter and downsample
    
    %Downsampling factors --? will be integral based on above computations
    ds_factor1 = (1000/isi1)/ DS_rate;
  
    %varargin: changing filter order and type for decimation
    %other options not yet added
    
   
    A1 = decimate( new_data1, ds_factor1, 'fir'); 
    T1 = decimate( T1, ds_factor1,'fir');

end