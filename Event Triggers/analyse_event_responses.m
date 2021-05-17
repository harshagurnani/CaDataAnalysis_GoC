function EventResposes = analyse_event_responses( Norm, Event, EventResponses, p, exp, Beh )
%% Isolate responses to a particular event occurring on all trials + corresponding behaviour

for sn = 1:exp.nSess
   sess = exp.session_types{exp.use_sessions(sn)};

   
   %% isolate matrix of puff responses and behaviour
   
   % empty for baseline sessions
   PuffResponses.(sess).n = [];
   PuffResponses.(sess).t = [];
   PuffResponses.(sess).AlignedN = [];
   PuffResponses.(sess).AlignedT = [];
   
   if ~isempty(Puff.(sess))
        
        extratime = [ -PuffResponses.params.Window(1)   PuffResponses.params.Window(2)];   %convert to TSOI extratime standard = [time befor event, time after event]
        
        % dF/F responses
        for roi=1:p.nROIs
            [temp.n{roi}, temp.t{roi}] = get_TSOI( Norm.(sess).n(:,roi), Norm.(sess).t(:,roi), Puff.(sess), extratime );
        end
        respsize = min(arrayfun( @(roi) length(temp.n{roi}{1}), 1:p.nROIs));        % min number of timepoints
        [PuffResponses.(sess).n,  PuffResponses.(sess).t , PuffResponses.(sess).pufftime] = deal(nan( respsize, p.nTrials, p.nROIs ));
        % crop responses of all ROI and trials to same number of timepoints
        for trial=1:p.nTrials
        for roi=1:p.nROIs
            PuffResponses.(sess).n(:,trial,roi)             = temp.n{roi}{trial}(1:respsize);                           % df/F in puff window 
            PuffResponses.(sess).t(:,trial,roi)             = temp.t{roi}{trial}(1:respsize);                           % actual timestamps of puff dF/F responses
            PuffResponses.(sess).pufftime(:, trial, roi)    = temp.t{roi}{trial}(1:respsize) - Puff.(sess)(trial,1);    % puff-aligned time
        end
        end
        
        % Behaviour - Speed
        [temp.loco, temp.locot] = get_TSOI( Speed.(sess)(:,2), Speed.(sess)(:,1), Puff.(sess), [500 1500]);
        respsize = min(arrayfun( @(trial) length(temp.locot{trial}), 1:p.nTrials));
        % crop speed of all trials to same number of timepoints
        for trial=1:p.nTrials
            PuffResponses.(sess).loco{trial}        = temp.loco{trial}(1:respsize);                                     % Encoder data around puff for each trial
            PuffResponses.(sess).locot{trial}       = temp.locot{trial}(1:respsize);                                    % timestamps of encoder data
            PuffResponses.(sess).loco_pufft{trial}  = temp.locot{trial}(1:respsize) - Puff.(sess)(trial,1);             % puff-aligned timestamps
        end
        
        % Behaviour - Whisking VidMI
        [temp.whisk, temp.whiskt] = get_TSOI( VidMI.(sess).(Whisk.vidseg)(:,1), VidMI.(sess).(Whisk.vidseg)(:,2), Puff.(sess), [500 1500]);
        respsize = min(arrayfun( @(trial) length(temp.whiskt{trial}), 1:p.nTrials));
        for trial = 1:p.nTrials
            PuffResponses.(sess).whisk{trial}       = temp.whisk{trial}(1:respsize);                                    % Whisker MI around puff for each trial
            PuffResponses.(sess).whiskt{trial}      = temp.whiskt{trial}(1:respsize);                                   % timestamps for whisker MI
            PuffResponses.(sess).whisk_pufft{trial} = temp.whiskt{trial}(1:respsize) - Puff.(sess)(trial,1);            % Puff-aligned timestamps

        end
        clear temp respsize
   
   
    %% Averaged responses to airpuff
    
        % Across all trials
        dt = 1000/PuffResponses.params.dsRate;
        
        locodt  = nanmean(diff(PuffResponses.(sess).locot{1}));
        whiskdt = nanmean(diff(PuffResponses.(sess).whiskt{1}));
        smoothbins.loco     = ceil( PuffResponses.params.smoothLoco/ locodt );
        smoothbins.whisk    = ceil( PuffResponses.params.smoothWhisk/whiskdt );
        
        mint = max(max(PuffResponses.(sess).pufftime(1, :, :))); maxt = min(min(PuffResponses.(sess).pufftime(end, :, :))); % shared interval = [mint, maxt]
        newint = [ceil(mint/dt), floor(maxt/dt)]; nt = newint(2)-newint(1) + 1; % number of timepoints
        newint = newint *dt;    %ms                                            % new interval - to ensure all ROI are covered
        
        
        PuffResponses.(sess).AlignedN                                           = nan(nt, p.nTrials, p.nROIs) ;
        [ PuffResponses.(sess).AlignedLoco, PuffResponses.(sess).AlignedWhisk ] = deal( nan(nt, p.nTrials) );
        
        for trial=1:p.nTrials
            for roi = 1:p.nROIs
                [ PuffResponses.(sess).AlignedN(:, trial, roi), PuffResponses.(sess).AlignedT ] = resample_ds( PuffResponses.(sess).n(:,trial,roi), PuffResponses.(sess).pufftime(:, trial, roi), ...
                                                                                                               PuffResponses.params.dsRate, newint );
                
                PuffResponses.(sess).baseline(trial, roi) = 0;
                b1 = find(PuffResponses.(sess).AlignedT>PuffResponses.params.baselineWindow(1),1); 
                if ~isempty(b1) 
                    b2 = find(PuffResponses.(sess).AlignedT<PuffResponses.params.baselineWindow(2));
                    if ~isempty(b2)
                        b2 = min(length(PuffResponses.(sess).AlignedT), b2(end)+1 );
                        PuffResponses.(sess).baseline(trial,roi) = nanmean( PuffResponses.(sess).AlignedN(b1:b2, trial, roi) );
                    end
                end
                PuffResponses.(sess).AlignedN(:, trial, roi) = PuffResponses.(sess).AlignedN(:, trial, roi) - PuffResponses.(sess).baseline(trial,roi);
                        
            end
            [ PuffResponses.(sess).AlignedLoco(:, trial), ~ ] = resample_ds( smooth(PuffResponses.(sess).loco{trial}, smoothbins.loco), PuffResponses.(sess).loco_pufft{trial}, ...
                                                                             PuffResponses.params.dsRate, newint );
            [ PuffResponses.(sess).AlignedWhisk(:, trial), ~ ] = resample_ds( smooth(PuffResponses.(sess).whisk{trial}, smoothbins.whisk), PuffResponses.(sess).whisk_pufft{trial}, ...
                                                                             PuffResponses.params.dsRate, newint );
        end
        
        % Across no locomotion trial
        PuffResponses.(sess).noLoco = nanmean(abs(PuffResponses.(sess).AlignedLoco)) < PuffResponses.params.noLoco_MeanE;
        
            
    
    end
   
end
% save([exp.analysed, '\PuffResponse.mat'],'PuffResponses','-v7.3')

end