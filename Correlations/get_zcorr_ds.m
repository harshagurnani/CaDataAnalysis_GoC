function [cc, rho] = get_zcorr_ds( data1, TS1, data2, TS2, DS_rate,varargin)
%% SYNTAX: [ cc, lags ] = get_zcorr_ds( data1, SR1, data2, SR2, DS_rate) OR
%                       
% ARGUMENTS:
%   - data1   [n x 1 (FLOAT)]:    
%           First data vector. Could be dFF  or speed or whisker position
%   - TS1     [n x 1 (FLOAT)]:       
%           Time stamps for first vector (ms)
%   - data2   [n x 1 (FLOAT)]:    
%           Second data vector. Could be dFF  or speed or whisker position.
%   - TS2     [n x 1 (FLOAT)]:
%           Time stamps of second vector (ms)
%   - DS_rate [INT]:
%           Rate of downsampled data for both vectors. (Hz) 
%
%   optional:
%
%
%   OUTPUT:
%
%       -
%
%
% Harsha Gurnani. March 2017

    %% same intervals for both data vectors
    same_interval=true;

    %% Interpolate data1 and data2
    %isi:   Mean 1000/Sample_Rate(Hz), ceiling to get new sampling
    isi1 = ceil( mean(TS1(2:end)-TS1(1:end-1) ));
    isi2 = ceil( mean(TS2(2:end) - TS2(1:end-1) ));
    %%%%%%%%% Ceil or floor?
    
    %fownsampling factors need to be integer so resample at one of the
    %divisor ISIs
    isi_to_use = divisors(1000/DS_rate);
    if (isi1 > 2*isi_to_use(end) ) || ( isi2 > 2*isi_to_use(end) )
        error('Downsampling rate is too high')
    end
    
    %Closest isi to use for resampling
    new_isi = find( isi_to_use>=isi1,1);
    if isempty(new_isi)
        isi1 = isi_to_use(end);
    else
        isi1=isi_to_use(new_isi);
    end
    new_isi = find( isi_to_use>=isi2,1);
    if isempty(new_isi)
        isi2 = isi_to_use(end);
    else
        isi2=isi_to_use(new_isi);
    end
    

    % Time interval: First and last time stamp
    if same_interval
        interval1 = [ max(TS1(1), TS2(1)), min(TS1(end),TS2(end)) ];
        interval2 = interval1;
    else
        interval1 = [TS1(1), TS1(end)];
        interval2 = [TS2(1), TS2(end)];
    end
    
    %Resample with spline interpolation
    interp_method = 'linear';
%     if any( mod(data1(2:end)-data1(1:end-1), 1 ) ~= 0 )
%         interp_method = 'spline';
%     end
    new_data1 = interp1( TS1, data1,  [interval1(1):isi1:interval1(2)] , interp_method);     % Linear is better if you don't want to introduce negative values in something that is supposed to be positive eg. probabilities becau
    
    interp_method = 'linear';
%     if any( mod(data2(2:end)-data2(1:end-1), 1 ) ~= 0 ),        interp_method = 'spline';    end
    new_data2 = interp1( TS2, data2,  [interval2(1):isi2:interval2(2)] , interp_method);

    
    %% Filter and downsample
    
    %Downsampling factors
    ds_factor1 = (1000/isi1)/ DS_rate;
    ds_factor2 = (1000/isi2)/ DS_rate;
    
    %varargin: changing filter order and type for decimation
    switch nargin
        case 5 %(default)
            A1 = decimate( new_data1, ds_factor1);
            A2 = decimate( new_data2, ds_factor2);
    end

    
    n1 = max(size(A1));
    n2 = max(size(A2));
    %reshaping to ensure they are vectors
    [cc, rho] = corrcoef(reshape(A1,[n1,1]), reshape(A2,[n2,1]));
    cc = cc(2);
    rho = rho(2);
end