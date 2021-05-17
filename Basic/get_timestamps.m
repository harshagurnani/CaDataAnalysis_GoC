function [data_time, params] = get_timestamps(params, ROIs_cycletime)
%% Returns timestamp matrix - same size as rawdata using given times for each point in cycle.
% Arguments: 
%   - params            Structure with fields like ROI_type, which ROIs/trials (not thoroughly tested for arrays of random ROIs but logically seems sound). 
%                       Other details read from params are size of data matrix - no. of timepoints, trials, ROIs
%   - ROIs_cycletime    2-column array. First column are acquisition times of successive points/lines in the cycle. Column 2 has just 1 value
%                       which is the total length of the cycle used to shift successive cycles.
%
% Returns:
%   - data_time         Timestamp matrix. 3D for pointscans and patch-averaged data. 5D for raw patch data.
%   - params            Updated params structure
%

    switch params.ROI_type
        case 'poi'
            %Same shape as raw_data
            data_time = nan(params.nTimepoints, params.nTrials, params.nROIs);

%             % single line acq time
%             t_per_line = params.nLinePixels*params.dwell_time;   % time to acquire last(single) line
%             %%% for new INtertrial FIFO, trial length can directly be calculated from trigger times - qqqqqqqq    
%             params.dwelltime_per_pixel_us = params.dwell_time; %0.0001
%             %% Check relative or HW timestamps
%             if exist([params.exp_path '/Single cycle relative times_HW.txt'], 'file')
%                 rel_time_given = false;
%                 line_times = importdata([params.exp_path '/Single cycle relative times_HW.txt']);
%                 if isstruct(line_times), line_times = line_times.data; end
%         %         exp_size = (params.maxLine * params.nTimepoints +1)* params.maxTrial; %expected no. of lines, with zero at reset of every trial
%                 % Remove rest zeros, so no +1 line per trial       
% %                 exp_size = (params.maxLine * params.nTimepoints)* params.maxTrial;
%                 %         line_times = line_times(end+1-exp_size:end,1:2)/1e3;     % timestamp for each line in ms         %to allow for initial zeros
%                 % Or With MC - tentative - qqqqq
%                 begin_time  = find(line_times(:,1)>0 ,1)-1;
%                 all_lines   = line_times(begin_time:end,:);
%                 MC_lines    = (all_lines(:,2)~=0);    %has non-zero time in col 2 (at reset is (0,0))
%                 line_times  = all_lines(all_lines(:,1)~=0,1)/1e3;    % ms (trial reset + roi line scans)
%                 MC_times    = all_lines( MC_lines,2)/1e3;    % ms (MC relative times per trial ) 
% %                 if length(line_times)~=exp_size, error('Did not find expected number of scan times '); end
%                 nTrialLines = length(line_times)/params.maxTrial;
%                 total_prev_im_time = [ 0; line_times((1:params.maxTrial-1)*nTrialLines)+t_per_line ];  %trial length = last imaging time + line acqtime! - qqq
%                 total_prev_im_time = cumsum(total_prev_im_time);                            %total time imaging before trials 1 to last (does not include intertrial time)
%             else,     rel_time_given = true;
%             end
            %Acquisition rate
            cycle_length = (ROIs_cycletime(1,2))/1000; % in ms.
            params.acquisition_rate = 1e3/cycle_length; % in Hz.

            %Multiple repetitions before movement correction
            params.Reps_Per_MC_Cycle = size(ROIs_cycletime,1)/params.maxROI;
            cycles_elapsed = floor( (0:params.nTimepoints-1)/params.Reps_Per_MC_Cycle );

            % Average acquisition rate - UPDATED for multiple acquisition
            % withing single MC cycle
            params.acquisition_rate = (1e3/cycle_length)*params.Reps_Per_MC_Cycle; % in Hz.

            %Block of relative cycle time for all ROIs and rep within cycle
            ROIs_reltimes = nan(params.Reps_Per_MC_Cycle, params.nROIs); 
            for rep = 1: params.Reps_Per_MC_Cycle
               ROIs_reltimes(rep,:) = reshape(ROIs_cycletime((rep-1)*params.maxROI + params.ROIs',1),[1,params.nROIs])/1000; %in ms 
            end

            %%Time between trials
            %--------- New update - checking initial junk by putting --------%
            %--------- limits of 2*trial length. Needs to be tested  --------%
            trial_length = params.nTimepoints/params.acquisition_rate;  %in second
            if params.add_intertime
                %-------- Comment this out if it's buggy:
                Time_between_trials = get_time_bw_trials([params.exp_path, '/Speed_Data/']);%, 2*trial_length);  %ms
                
                %-------- Uncomment this to reuse older, working version:
%               Time_between_trials = get_time_bw_trials([params.exp_path, '/Speed_Data/']);  %ms
                if isempty(Time_between_trials)
                    warning('Time between trials is empty, dividing dead space equally')
                    sd = dir([params.exp_path, '/Speed_Data/Speed data*.txt']); speed_data = importdata([ sd.folder, '/', sd.name]);
                    clock_time      = (speed_data.data(end,1)-speed_data.data(1,1))*0.001;   %ms
                    imaging_time    = trial_length*1000 *params.nTrials*1000;                %ms
                    Time_between_trials = repmat( (clock_time-imaging_time)/(params.nTrials-1), 1, params.nTrials-1);
                end
    
            end

            for nt = 1:params.nTrials
                trial_n = params.trials(nt);
                total_prevtime = (trial_n-1)*params.nTimepoints*cycle_length ;  %Total elapsed time (in ms)
                if trial_n > 1 && exist('Time_between_trials', 'var')
                    total_prevtime = total_prevtime + sum(Time_between_trials(1:trial_n-1));
                end
                data_time(:,nt,:) = reshape( ...
                  repmat( cycle_length *cycles_elapsed', [1,params.nROIs]) + ... cycles elapsed * cycle_length
                  repmat( ROIs_reltimes, [params.nTimepoints/params.Reps_Per_MC_Cycle, 1] ) + ... relative time of POI per cycle
                  repmat(total_prevtime, [params.nTimepoints, params.nROIs]), ... elapsed time
                  [params.nTimepoints,1,params.nROIs]);
            end
            
            %%
        case 'patch'
            if exist([params.exp_path '/Single cycle relative times_HW.txt'], 'file'),
                line_times = importdata([params.exp_path '/Single cycle relative times_HW.txt']);
                if isstruct(line_times), line_times = line_times.data; end
                exp_size = params.nPatchLines * params.nROIs * params.nTimepoints * params.nTrials;
                line_times = line_times(end+1-exp_size:end,1)/1000;     % timestamp for each line in ms
                temp_speed = get_speed_data([params.exp_path, '/Speed_Data/']);
                total_prev_time = [0 arrayfun( @(nt) temp_speed{nt}(end,1),1:params.maxTrials )];
                clear temp_speed
            end
            
            if params.average_patch
             data_time = nan(params.nTimepoints, params.nTrials, params.nROIs );
            else
             data_time = nan(params.nLinePixels, params.nPatchLines, params.nTimepoints, params.nTrials, params.nROIs );
            end

            nLines = params.nPatchLines*params.nROIs;
            
            %Acquisition rate
            if length(ROIs_cycletime)<nLines
                cycle_length = line_times(nLines+1)-line_times(1);    
                rel_time_given  = false;
            else
                cycle_length = (ROIs_cycletime(1,2))/1000; % in ms.
                rel_time_given = true;
            end
            params.acquisition_rate = 1e3/cycle_length; % in Hz.
            
            %Fixed rates (check with header file):
            params.dwelltime_per_pixel_us = 0.100;    %in us
            
            t_per_line = params.nLinePixels * params.dwelltime_per_pixel_us;
            
            %%%%%%%%%%%%%%% To add reps per MC for patches - Usually 1
            %%%%%%%%%%%%%%% unless cyclength is less than 2 ms.
            %params.Reps_Per_MC_Cycle = size(ROIs_cycletime,1)/(params.maxROI*params.nPatchLines);
            %%%%%%%%%%%%%%%
            cycles_elapsed = floor( (0:params.nTimepoints-1) );
            %%Time between trials
            trial_length = params.nTimepoints/params.acquisition_rate;  %in second

            if params.add_intertime
                Time_between_trials = get_time_bw_trials([params.exp_path, '/Speed_Data/']);
                if isempty(Time_between_trials)
                    % in case there were no intertrial triggers and the
                    % encoder data was recorded continuously. Last
                    % timestamp = total recording time.
                    warning('Time between trials is empty, dividing dead space equally')
                    sd = dir([params.exp_path, '/Speed_Data/Speed data*.txt']); speed_data = importdata([ sd.folder, '/', sd.name]);
                    clock_time      = (speed_data.data(end,1)-speed_data.data(1,1))*0.001;   %ms
                    imaging_time    = trial_length*1000 *params.nTrials;                %ms
                    Time_between_trials = repmat( (clock_time-imaging_time)/(params.nTrials-1), 1, params.nTrials-1);
                end

            end

            if rel_time_given
                %Block of relative cycle time for all ROIs and rep within cycle
                ROIs_reltimes = nan(params.nLinePixels, params.nPatchLines, params.nROIs); 
                for line_n =1:params.nPatchLines
                    for roi_n = 1:params.nROIs
                        ROIs_reltimes(:,line_n,roi_n) = reshape( (0:params.nLinePixels-1) * params.dwelltime_per_pixel_us  + ... Pixel number*dwell time along line
                                                            ROIs_cycletime(line_n+(roi_n-1)*params.nPatchLines, 1), ... % Starting of line  
                                                        [params.nLinePixels, 1, 1] ) /1000; %in ms 

                    end
                end


                inter_trial_time = 0;       % Changed later if asked to add inter-trial time
                for nt = 1:params.nTrials
                    trial_n = params.trials(nt);
                    total_prevtime = (trial_n-1)*params.nTimepoints*cycle_length ;  %Total elapsed time during IMAGING ONLY (in ms)
                    %%% previous code was fine, but i thought i was updating
                    %%% prev_time outside of loop  and it would fail for trials
                    %%% listed in random (ie. non-monotonic order - hence
                    %%% separated to prevent future panic)
                    if trial_n>1 && exist('Time_between_trials', 'var')
                        inter_trial_time = sum(Time_between_trials(1:trial_n-1));
                    end

                    if params.average_patch
                        reltime = reshape( squeeze( mean( mean( ROIs_reltimes ) )  ), ...
                                            [1, params.nROIs] );


                        data_time(:,nt,:) = reshape( ...
                          repmat( cycle_length *cycles_elapsed', [1,params.nROIs]) +                   ... cycles elapsed * cycle_length
                          repmat( reltime,                       [params.nTimepoints, 1]) +            ... relative time of POI per cycle
                          repmat( total_prevtime,                [params.nTimepoints, params.nROIs]) + ... elapsed imaging time
                          repmat( inter_trial_time,              [params.nTimepoints, params.nROIs]),  ... inter-trial dead time
                          [params.nTimepoints,1,params.nROIs]);
                    else
                        reltime = reshape( ROIs_reltimes, [params.nLinePixels, params.nPatchLines, 1, 1, params.nROIs]);


                        data_time(:,:,:,nt,:) = reshape( ...
                          repmat( reshape(cycle_length *cycles_elapsed', [1,1,params.nTimepoints]),...
                                  [params.nLinePixels, params.nPatchLines, 1, 1, params.nROIs]) + ... cycles elapsed * cycle_length
                          repmat( reltime, [1,1,params.nTimepoints, 1, 1] ) + ... relative time of POI per cycle
                          total_prevtime, ... elapsed time
                          [params.nLinePixels, params.nPatchLines, params.nTimepoints,1,params.nROIs]);
                    end
                end
            
            else
                nTrialLines = length(line_times)/params.nTrials;
                    
                for nt = 1:params.nTrials
                    trial_n = params.trials(nt);
                    trialtime = line_times((trial_n-1)*nTrialLines+1:trial_n*nTrialLines);
                    inter_trial_time = 0;       % Changed later if asked to add inter-trial time
                    if trial_n>1 && exist('Time_between_trials', 'var')
                        inter_trial_time = Time_between_trials(trial_n-1);
                    end
                    if params.average_patch                  
                        trialtime = reshape( trialtime, [params.nPatchLines, params.nROIs, params.nTimepoints] ); 
                        trialtime = squeeze(nanmean(trialtime,1)) ;
                        trialtime = permute(trialtime, [2 1]);
                        data_time(:, nt, :) =  trialtime+ ...
                                               repmat( total_prev_time(trial_n) + inter_trial_time, size(trialtime)); %ms
                    else
                        trialtime = reshape( trialtime, [1, params.nPatchLines, 1, params.nROIs, params.nTimepoints] );
                        trialtime = repmat( trialtime, [params.nLinePixels, 1, 1, 1, 1] );
                        trialtime = permute(trialtime, [1 2 5 3 4]);
                        data_time(:, :, :, nt, :) = trialtime + ...
                                               repmat( total_prev_time(trial_n)+ inter_trial_time, size(trialtime)); %ms
                    end
                    
                end
                                
            end
    end
        
end