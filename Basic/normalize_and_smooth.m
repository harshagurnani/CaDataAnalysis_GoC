function [norm_data, baseline] = normalize_and_smooth(rawdata, timedata, params, varargin)
%% Normalize and smooth: calculate dF/F and smoothen. Also returns baseline used for normalisation for each ROI(trial)
% Arguments:
%   - rawdata       3D array (Time,Trial,ROI) of raw F
%   - timedata      Acquisition Timestamps(ms) 
%                   Same size as raw data. Used to determine bins for
%                   smoothing based on smoothing window specified in (ms)
%   - params        Structure with normalisation parameters in fields
%                   'norm_method', 'baseline_per_trial', 'base_percentile',
%                   smooth_type', 'smooth_scale'. 
%                   If it is empty, then default values are used.
% 
% SMOOTHING WINDOW
%   -regular - BOOL argument - for equal time steps between consecutive data points (except trial changes, which is automatically adjusted for.) 
%              True by default.
%
% Returns:
%   - norm_data     Normalised (and smoothed) data (dF/F) - Same size as raw data.
%   - baseline      Array (ROIxTrial) 
%
%       IMPORTANT NOTE: All temporal smoothing is now causal i.e. only
%       measurements in smoothing window BEFORE current timepoint are
%       averaged over - future measurements are not impacted. Important in
%       cases of latency calculations to avoid spillover due to smoothing.
%--------------------------


    % regular = true implies data at regular, constant intervals. 
    if nargin == 3
        regular = true;
    elseif (nargin == 4 && islogical(varargin{1}) )
        regular = varargin{1};
    else
        error('Only accept 3 (rawdata, timedata, params) or 4 arguments (optional argument must be logical)')
    end
    
    
    tp = size(rawdata,1);   %For eg. with concatenated or smoothed data
    poi = size(rawdata,3);
    trial = size(rawdata,2);

    
    %%%% when params is empty, add default values.
    if isempty( params)
        disp('No normalisation /smoothing params given. Using default values.')
        params.norm_method          = 'percentile';
        params.base_percentile      = 10;
        params.baseline_per_trial   = false;
        params.smooth_type          = 'Moving average';
        params.smooth_scale         = 250;
        try
            params.acquisition_rate     = timedata(2,1)-timedata(1,1);
        catch
            params.acquisition_rate     = 50;
            warning('No acquisition rate known - smoothing over 5 time bins' )
        end
        if trial > 1
            warning('More than 1 trial in array - assuming they are continuous' )
            params.trials = 1:trial;
        end
    end
    
    if trial>1 && isfield(params, 'trials')
        diffs = params.trials(2:end)-params.trials(1:end-1);
        continuous = ~any( diffs ~=1 ) || ~any( diffs ~=-1 );
    else
        continuous = true;
    end
   
    %%% causal filter for smoothing?
    if ~isfield( params, 'smooth_causal')
        params.smooth_causal = false;
    end
    
    %-------------------------------------------------------------------%
    %% Normalization
    %-------------------------------------------------------------------%
    if params.baseline_per_trial
        switch params.norm_method
            case 'median'
                baseline = median(rawdata,1);
            case 'percentile'
                baseline = prctile(rawdata, params.base_percentile, 1);
        end
        F0 = repmat(baseline, tp, 1);
        norm_data = (rawdata - F0)./F0;
        
    else
        baseline = nan(1,poi);
        parfor pn=1:poi
             temp=reshape(rawdata(:,:,pn), [tp*trial,1]); 
             %%%Add exp-fit fto correct for bleaching
             switch params.norm_method
                 case 'median'
                    baseline(pn) = median(temp);
                 case 'percentile'
                     baseline(pn) = prctile(temp, params.base_percentile);
             end
        end
        F0 = repmat(reshape(baseline,[1,1,poi]), [tp,trial,1]);
        norm_data = (rawdata - F0)./F0;

    end
    
    %-------------------------------------------------------------------%
    %% Smoothening normalized data
    %-------------------------------------------------------------------%
    smooth_bins = ceil(params.acquisition_rate * params.smooth_scale / 2000);
    if smooth_bins > 0

        if regular && continuous
        % Continuous trials, and all intratrial points are regularly spaced.
 
%------------------- OLD, INEFFICIENT CODE -- to be deleted later ---------------------------%        
            % Pad inter-trial edges for smoothing.
%             if trial > 1
%                 temp1 = cat( 1, repmat(norm_data(1, 1, :), [smooth_bins, 1, 1] ),...
%                                     norm_data(:,1,:), ...
%                                    (norm_data(tp,1,:)+norm_data(1,2,:))/2,...
%                                     norm_data(1:smooth_bins-1, 2, :) );
%                 temp_mid = nan(tp+2*smooth_bins,trial-2,poi);
%                 for tn=2:trial-1
%                     temp_mid(:,tn-1,:) = cat(  1, norm_data(tp-smooth_bins+2:tp, tn-1, :), ...
%                                             (norm_data(tp,tn-1,:) + norm_data(1,tn,:))/2, ...
%                                              norm_data(:,tn,:), ...
%                                             (norm_data(tp,tn,:) + norm_data(1,tn+1,:))/2, ...
%                                              norm_data(1:smooth_bins-1, tn+1, :) );
% 
%                 end
%                 temp_end = cat(1, norm_data(tp-smooth_bins+2:tp,trial-1,:), ...
%                                (norm_data(tp,trial-1,:)+norm_data(1,trial,:))/2,...
%                                 norm_data(:, trial, :) ,...
%                                 repmat(norm_data(tp, trial, :), [smooth_bins, 1, 1] ));
%                 temp = cat(2, temp1, temp_mid, temp_end);
%             else
%                 temp = cat( 1, repmat(norm_data(1, 1, :), [smooth_bins, 1, 1] ),...
%                                    norm_data(:,1,:),...
%                             repmat(norm_data(tp, 1, :), [smooth_bins, 1, 1] ) );
%             end
% 
%             switch params.smooth_type
%             	case 'Moving average'
%                     for tn=1:trial
%                     parfor pn=1:poi
%                         expanded = smooth( temp(:,tn,pn), 2*smooth_bins+1 );
%                         norm_data(:,tn,pn) = expanded(smooth_bins+1:tp+smooth_bins);
%                     end
%                     end
%                 case 'Moving median'
%                     for tn=1:trial
%                     parfor pn=1:poi
%                     	expanded = medfilt1( temp(:,tn,pn), 2*smooth_bins+1 );
%                         norm_data(:,tn,pn) = expanded(smooth_bins+1:tp+smooth_bins);
%                     end
%                     end
%             end
%       
%-------------------------------------------------------------------------%

            % Continuous trials, and all intratrial points are regularly spaced.
        
            temp = reshape( norm_data, [tp*trial, poi] );
            if params.smooth_causal
                % shift data by smooth_bins
                temp = [repmat( temp(1,:), [smooth_bins, 1]);
                        temp ];
            end
            
            switch params.smooth_type
            	case 'Moving average'
                    parfor pn=1:poi
                        temp(:,pn) = smooth( temp(:,pn), 2*smooth_bins+1 );                        
                    end
                    
                    
                case 'Moving median'
                    
                    parfor pn=1:poi
                    	temp(:,pn) = medfilt1( temp(:,pn), 2*smooth_bins+1 );                                           
                    end
            end
            if params.smooth_causal
                % Remove last smooth_bins timepoints
                temp = temp(1:end-smooth_bins, :);
            end
            norm_data = reshape( temp, [tp, trial, poi] ); 
            clear temp


        elseif regular  
            %all intratrial points are regularly spaced.
            % But with missing/unsorted trials - each trial is individually
            % padded and smoothed, thus all selected trials will have noisy ends
            
            switch params.smooth_type
            	case 'Moving average'
                    for tn=1:trial
                    for pn=1:poi
                         temp = smooth( [repmat(norm_data(1,tn,pn), [smooth_bins,1,1]);...
                                                norm_data(:,tn,pn)],          2*smooth_bins+1);
                         norm_data(:,tn,pn) = temp(1:end-smooth_bins);
                    end
                    end
                case 'Moving median'
                    for tn=1:trial
                    for pn=1:poi
                         temp = medfilt1( [repmat( norm_data(1,tn,pn), [smooth_bins,1,1]);...
                                                   norm_data(:,tn,pn)],          2*smooth_bins+1);
                         norm_data(:,tn,pn) = temp(1:end-smooth_bins);                        
                    end
                    end
            end

        else
            % timepoints/trials are not regularly spaced or sorted - each
            % timepoint is individually smoothed with measurements within a
            % small time window before it.
            
            temp=norm_data;
            switch lower( params.smooth_type )
            	case 'moving average'
                    for tn=1:trial
                    parfor pn=1:poi
                    for dt = 1:tp
                        tr = ( ( timedata(dt,tn,pn) - timedata(:,tn,pn)  <= params.smooth_scale) &...
                               ( timedata(:,tn,pn) -  timedata(dt,tn,pn) <=0  )                         );
                        norm_data(dt,tn,pn) = mean( temp(tr,tn,pn) );
                    end
                    end
                    end
                case 'moving median'
                    for tn=1:trial
                    parfor pn=1:poi
                    for dt = 1:tp
                        tr = ( ( timedata(dt,tn,pn) - timedata(:,tn,pn)  <= params.smooth_scale) &...
                               ( timedata(:,tn,pn) -  timedata(dt,tn,pn) <=0  )                        );
                        norm_data(dt,tn,pn) = median( temp(tr,tn,pn) );
                    end
                    end
                    end
            end
            clear temp
            
            
        end
    
    end

            
end
