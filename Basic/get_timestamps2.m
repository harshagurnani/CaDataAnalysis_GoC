function [data_time, params] = get_timestamps2(params, ROIs_cycletime)
%% Returns timestamp matrix - same size as rawdata
% Arguments: 
%   - params            Structure with fields like ROI_type, which ROIs/trials (not thoroughly tested for arrays of random ROIs but logically seems sound). 
%                       Other details read from params are size of data matrix - no. of timepoints, trials, ROIs,
%                       Also encoder_shift_dt (back compatibility for missing encoder time)
%   - ROIs_cycletime    2-column array. First column are acquisition times of successive points/lines in the cycle. Column 2 has just 1 value
%                       which is the total length of the cycle used to shift successive cycles.
%                       From Single cycle relative times.txt (back compatibility)
%
% Returns:
%   - data_time         Timestamp matrix. 3D for pointscans and patch-averaged data. 5D for raw patch data.
%   - params            Updated params structure
% HG, March 2019

    % single line acq time
    % fill in missing fields - qqq temporary
    params.maxLine = params.maxROI*params.nPatchLines;
    params.maxTrial = params.maxTrials;
    params.encoder_shift_dt = 0;
    
    t_per_line = params.nLinePixels*params.dwell_time;   % time to acquire last(single) line
    %%% for new INtertrial FIFO, trial length can directly be calculated from trigger times - qqqqqqqq    
    params.dwelltime_per_pixel_ms = params.dwell_time; %0.0001
    %% Check relative or HW timestamps
    if exist([params.exp_path '/Single cycle relative times_HW.txt'], 'file')
        rel_time_given = false;
        line_times = importdata([params.exp_path '/Single cycle relative times_HW.txt']);
        if isstruct(line_times), line_times = line_times.data; end
%         exp_size = (params.maxLine * params.nTimepoints +1)* params.maxTrial; %expected no. of lines, with zero at reset of every trial
        % Remove rest zeros, so no +1 line per trial       
        exp_size = (params.maxLine * params.nTimepoints)* params.maxTrial;
        
        %         line_times = line_times(end+1-exp_size:end,1:2)/1e3;     % timestamp for each line in ms         %to allow for initial zeros
        % Or With MC - tentative - qqqqq
        begin_time  = find(line_times(:,1)>0 ,1)-1;
        all_lines   = line_times(begin_time:end,:);
        MC_lines    = (all_lines(:,2)~=0);    %has non-zero time in col 2 (at reset is (0,0))
        line_times  = all_lines(all_lines(:,1)~=0,1)/1e3;    % ms (trial reset + roi line scans)
        MC_times    = all_lines( MC_lines,2)/1e3;    % ms (MC relative times per trial ) 
        if length(line_times)~=exp_size, error('Did not find expected number of scan times '); end
        nTrialLines = length(line_times)/params.maxTrial;
        total_prev_im_time = [ 0; line_times((1:params.maxTrial-1)*nTrialLines)+t_per_line ];  %trial length = last imaging time + line acqtime! - qqq
        total_prev_im_time = cumsum(total_prev_im_time);                            %total time imaging before trials 1 to last (does not include intertrial time)
    else,     rel_time_given = true;
    end
    
    %% Acquisition rate
    if length(ROIs_cycletime)<params.maxLine
        % compute mean cycle length from first cycle of each trial
        cycle_length = arrayfun( @(jj) line_times(1+params.maxLine+(params.maxLine*params.nTimepoints+1)*(jj-1)),1:params.maxTrial );
        cycle_length = nanmean(cycle_length);   
%         rel_time_given  = false; %%% give priority to hw existing, not single cycle
        params.acquisition_rate = 1e3/cycle_length; % in Hz.

    else
        cycle_length = (ROIs_cycletime(1,2))/1000; % in ms.
%         rel_time_given = true;

        %Multiple repetitions before movement correction
        params.Reps_Per_MC_Cycle = size(ROIs_cycletime,1)/params.maxLine;
        cycles_elapsed = floor( (0:params.nTimepoints-1)/params.Reps_Per_MC_Cycle );
        params.acquisition_rate = (1e3/cycle_length)*params.Reps_Per_MC_Cycle; % in Hz.

    end

    
    %% Inter-trial time
    trial_length = params.nTimepoints/params.acquisition_rate * 1e3;  %in ms
    inter_trial_time = zeros( params.maxTrial, 1);  %zero by default
                    
    if params.add_intertime
        [~,Time_between_trials] = get_intertrial_time2(params);
        if isempty(Time_between_trials)
            % Bug check: in case there were no intertrial triggers and the  encoder data was recorded continuously. 
            warning('Time between trials is empty, dividing dead space equally')
            sd = dir([params.exp_path, '/Speed_Data/Speed data*.txt']); speed_data = importdata([ sd.folder, '/', sd.name]);
            clock_time      = (speed_data.data(end,1)-speed_data.data(1,1))*1e-3;   %ms    % Last timestamp = total recording time.
            imaging_time    = trial_length *params.maxTrial;                       %ms
            Time_between_trials = repmat( (clock_time-imaging_time)/(params.maxTrial-1), params.maxTrial-1, 1);
            Time_between_trials = Time_between_trials + params.encoder_shift_dt;%add missing time from encoder reset (back-compatibility)
        
        end
        inter_trial_time = cumsum( [0; Time_between_trials] );              %total intertrial non-imaging time

    end

    %% Imaging mode
    switch lower(params.ROI_type)            
        case 'patch'
            if params.average_patch
             data_time = nan(params.nTimepoints, params.nTrials, params.nROIs );
            else
             data_time = nan(params.nLinePixels, params.nPatchLines, params.nTimepoints, params.nTrials, params.nROIs );
            end

            if ~rel_time_given
            %% Use FIFO (hardware timestamps ) -> check for MC            
                % without MC
%                 if ~ params.online_MC
                    
                for nt = 1:params.nTrials
                    trial_n = params.trials(nt);        % to allow for non-cont trial selection
%                     trialtime = line_times((trial_n-1)*nTrialLines+2:trial_n*nTrialLines);
                    % removed the reset zero line
                    trialtime = line_times((trial_n-1)*nTrialLines+1:trial_n*nTrialLines);
                    if params.average_patch                  
                        trialtime = reshape( trialtime, [params.nPatchLines, params.nROIs, params.nTimepoints] ); 
                        trialtime = squeeze(nanmean(trialtime,1)) ;
                        trialtime = permute(trialtime, [2 1]);
                        data_time(:, nt, :) =  trialtime+ ...
                                               repmat( total_prev_im_time(trial_n) + inter_trial_time(trial_n), size(trialtime)); %ms
                    else
                        trialtime = reshape( trialtime, [1, params.nPatchLines, 1, params.nROIs, params.nTimepoints] );
                        trialtime = repmat( trialtime, [params.nLinePixels, 1, 1, 1, 1] );
                        trialtime = permute(trialtime, [1 2 5 3 4]);
                        data_time(:, :, :, nt, :) = trialtime + ...
                                                    repmat( (0:params.nLinePixels-1)'*params.dwelltime_per_pixel_ms, [1, params.nPatchLines, params.nTimepoints, 1, params.nROIs]) + ... Time along linescan
                                                    repmat( total_prev_im_time(trial_n)+ inter_trial_time(trial_n), size(trialtime)); %ms   
                    end                    
                end
%                 
%                 else
%                     % FIFO times for MC not written - qqqqqq
%                 end
                
                
            else
            %% Old format - Use relative times per cycle to generate timestamps
            
                %Block of relative cycle time for all ROIs and rep within cycle
                ROIs_reltimes = nan(params.nLinePixels, params.nPatchLines, params.nROIs); 
                for line_n =1:params.nPatchLines
                    for roi_n = 1:params.nROIs
                        ROIs_reltimes(:,line_n,roi_n) = reshape( (0:params.nLinePixels-1) * params.dwelltime_per_pixel_ms  + ... Pixel number*dwell time along line
                                                            ROIs_cycletime(line_n+(roi_n-1)*params.nPatchLines, 1), ... % Starting of line  
                                                        [params.nLinePixels, 1, 1] ) /1e3; %in ms 
                    end
                end

                for nt = 1:params.nTrials
                    trial_n = params.trials(nt);
                    total_prev_imtime = (trial_n-1)*params.nTimepoints*cycle_length ;  %Total elapsed time during IMAGING ONLY (in ms)
                    if params.average_patch
                        %t, tr, roi
                        reltime = reshape( squeeze( mean( mean( ROIs_reltimes,2 ),1 )  ), ...
                                            [1, params.nROIs] );
                                        
                        data_time(:,nt,:) = reshape( ...
                          repmat( cycle_length *cycles_elapsed', [1,params.nROIs]) +                   ... cycles elapsed * cycle_length
                          repmat( reltime,                       [params.nTimepoints, 1]) +            ... relative time of POI per cycle
                          repmat( total_prev_imtime,             [params.nTimepoints, params.nROIs]) + ... elapsed imaging time
                          repmat( inter_trial_time(trial_n),     [params.nTimepoints, params.nROIs]),  ... inter-trial dead time
                          [params.nTimepoints,1,params.nROIs]);
                    else
                        % full x,y,t, tr, roi
                        reltime = reshape( ROIs_reltimes, [params.nLinePixels, params.nPatchLines, 1, 1, params.nROIs]);

                        data_time(:,:,:,nt,:) = reshape( ...
                          repmat( reshape(cycle_length *cycles_elapsed', [1,1,params.nTimepoints]),...
                                  [params.nLinePixels, params.nPatchLines, 1, 1, params.nROIs])  + ... cycles elapsed * cycle_length
                          repmat( reltime,                      [1,1,params.nTimepoints, 1, 1] ) + ... relative time of POI per cycle
                          total_prev_imtime                                                      + ... elapsed time
                          inter_trial_time(trial_n),                                               ... intertrial dead time
                          [params.nLinePixels, params.nPatchLines, params.nTimepoints,1,params.nROIs]);
                    end
                end
          
            end
            
        case 'poi'
            
            %Same shape as raw_data
            data_time = nan(params.nTimepoints, params.nTrials, params.nROIs);
           
            if rel_time_given
            %% Relative times per cycle
            
                %Block of relative cycle time for all ROIs and rep within cycle
                ROIs_reltimes = nan(params.Reps_Per_MC_Cycle, params.nROIs); 
                for rep = 1: params.Reps_Per_MC_Cycle
                   ROIs_reltimes(rep,:) = reshape(ROIs_cycletime((rep-1)*params.maxROI + params.ROIs',1),[1,params.nROIs])/1000; %in ms 
                end

                for nt = 1:params.nTrials
                    trial_n = params.trials(nt);
                    total_prev_imtime = (trial_n-1)*params.nTimepoints*cycle_length ;  %Total elapsed time (in ms)
                    data_time(:,nt,:) = reshape( ...
                      repmat( cycle_length *cycles_elapsed',                [1,params.nROIs]) + ... cycles elapsed * cycle_length
                      repmat( ROIs_reltimes,                                [params.nTimepoints/params.Reps_Per_MC_Cycle, 1] ) + ... relative time of POI per cycle
                      repmat(total_prev_imtime + inter_trial_time(trial_n), [params.nTimepoints, params.nROIs]), ... elapsed time
                      [params.nTimepoints,1,params.nROIs]);
                end
            
            else
            %% READ from FIFO - to be tested -qqqqqq
                if ~ params.online_MC

                    for nt = 1:params.nTrials
                        trial_n = params.trials(nt);
%                         trialtime = line_times((trial_n-1)*nTrialLines+2:trial_n*nTrialLines);
                        trialtime = line_times((trial_n-1)*nTrialLines+1:trial_n*nTrialLines);  %removed the trial reset zero line
                                    
                        trialtime = reshape( trialtime, [params.nROIs, 1, params.nTimepoints] ); 
                        trialtime = permute(trialtime, [3 2 1]);
                        data_time(:, nt, :) =  trialtime+ ...
                                               repmat( total_prev_im_time(trial_n) + ...
                                                       inter_trial_time(trial_n),   size(trialtime)); %ms
                    end
                
                else
                    % FIFO times for MC not written - qqqqqq
                end
                
            end

    end
        
end