function [EventResponses] = get_event_responses2( EventPOI, Traces, TracesTime, extratime , Beh, EventParams )

   %% isolate matrix of puff responses and behaviour

   % empty for baseline sessions
   EventResponses.n = [];
   EventResponses.t = [];
   EventResponses.AlignedN = [];
   EventResponses.AlignedT = [];
   nEvents = size(EventPOI,1);
   
   whisk_mi = false;

    dt = 1000/EventParams.dsRate;

    temp={};nROI = size( Traces, 2); if size(TracesTime,2) == 1, TracesTime = repmat(TracesTime,1, nROI); end
    % dF/F responses
    for roi=1:nROI
        [temp.n{roi}, temp.t{roi}] = get_TSOI( Traces(:,roi), TracesTime(:,roi), EventPOI, extratime );
    end
    
    respsize = max(arrayfun( @(roi) length(temp.n{roi}{1}), 1:nROI));        % max number of timepoints
    [EventResponses.n,  EventResponses.t , EventResponses.eventtime] = deal(nan( respsize, nEvents, nROI ));
    % crop responses of all ROI and trials to same number of timepoints
    
    for trial=1:nEvents
    for roi=1:nROI
        respsize = length(temp.n{roi}{trial});
        EventResponses.n(1:respsize,trial,roi)             = temp.n{roi}{trial};                           % df/F in puff window 
        EventResponses.t(1:respsize,trial,roi)             = temp.t{roi}{trial};                           % actual timestamps of puff dF/F responses
        EventResponses.eventtime(1:respsize, trial, roi)   = temp.t{roi}{trial} - EventPOI(trial,1);       % puff-aligned time
    end
    end

    % Behaviour
    
    if ~isempty(Beh)
       beh_types = fieldnames(Beh);
       for bb=1:length(beh_types)
        beh_id = beh_types{bb};
        behdt  = nanmean(diff(Beh.(beh_id)(:,1)));
        smoothbins.(beh_id)     = ceil( EventParams.smoothBeh(bb)/ behdt );

        [temp.(beh_id).x, temp.(beh_id).t] = get_TSOI( smooth(Beh.(beh_id)(:,2),smoothbins.(beh_id)), Beh.(beh_id)(:,1), EventPOI, extratime);
        respsize = (arrayfun( @(trial) length(temp.(beh_id).x{trial}), 1:nEvents)); respsize=max(respsize(respsize>10));
        % crop speed of all trials to same number of timepoints
        for trial=1:nEvents
            try
            nt = length(EventResponses.(beh_id).x{trial}); 
            EventResponses.(beh_id).x{trial} = [temp.(beh_id).x{trial}; nan(respsize-nt,1)];
            EventResponses.(beh_id).t{trial} = [temp.(beh_id).t{trial}; nan(respsize-nt,1)];                                     % Encoder data around puff for each trial
            EventResponses.(beh_id).event_t{trial} = EventResponses.(beh_id).t{trial} - EventPOI(trial,1);             % puff-aligned timestamps
            catch
            warning(sprintf('Trial %d did not have sufficient speed data!',trial))
            EventResponses.(beh_id).x{trial}       = temp.(beh_id).x{trial};                                    % Encoder data around puff for each trial
            EventResponses.(beh_id).t{trial}       = temp.(beh_id).t{trial};                                  % timestamps of encoder data
            EventResponses.(beh_id).event_t{trial} = EventResponses.(beh_id).t{trial} - EventPOI(trial,1);            % puff-aligned timestamps
            end
        end

          
       end
    end
    clear temp respsize


%% Averaged responses to airpuff
    
% Across all trials
    dt = 1000/EventParams.dsRate;
    
    mint = nanmin(nanmin(EventResponses.eventtime(1,:,:))); 
    maxt = nanmax(nanmax(EventResponses.eventtime(end,:,:))); % max interval = [mint, maxt]
    newint = [ceil(mint/dt), floor(maxt/dt)]; nt = newint(2)-newint(1) + 1; % number of timepoints
    newint = newint *dt;    %ms                                            % new interval - to ensure all ROI are covered


    EventResponses.AlignedN                                           = nan(nt, nEvents, nROI) ;
    if ~isempty(Beh)
       beh_types = fieldnames(Beh);
       for bb=1:length(beh_types)
        beh_id = beh_types{bb};
        [ EventResponses.(beh_id).AlignedX ] = nan(nt, nEvents);
       end
    end
    
    for trial=1:nEvents
        for roi = 1:nROI
            [ EventResponses.AlignedN(:, trial, roi), EventResponses.AlignedT ] = resample_ds( inpaint_nans(EventResponses.n(:,trial,roi)), inpaint_nans( EventResponses.eventtime(:, trial, roi)), ...
                                                                                                           EventParams.dsRate, newint );

            EventResponses.baseline(trial, roi) = 0;
            b1 = find(EventResponses.AlignedT>EventParams.baselineWindow(1),1); 
            if ~isempty(b1) 
                b2 = find(EventResponses.AlignedT<EventParams.baselineWindow(2));
                if ~isempty(b2)
                    b2 = nanmin(length(EventResponses.AlignedT), b2(end)+1 );
                    EventResponses.baseline(trial,roi) = nanmean( EventResponses.AlignedN(b1:b2, trial, roi) );
                end
            end
            EventResponses.AlignedN(:, trial, roi) = EventResponses.AlignedN(:, trial, roi) - EventResponses.baseline(trial,roi);

        end
        if ~isempty(Beh)
        beh_types = fieldnames(Beh);
        for bb=1:length(beh_types)
        beh_id = beh_types{bb};
            
        try
        [ EventResponses.(beh_id).AlignedX(:, trial), ~ ] = resample_ds( EventResponses.(beh_id).x{trial}, EventResponses.(beh_id).event_t{trial}, ...smooth(EventResponses.loco{trial}, smoothbins.loco)
                                                                         EventParams.dsRate, newint );
        catch
        end
        end
        end

    % Across no locomotion trial
%     EventResponses.noLoco = nanmean(abs(EventResponses.AlignedLoco)) < EventResponses.params.noLoco_MeanE;


    end
end
