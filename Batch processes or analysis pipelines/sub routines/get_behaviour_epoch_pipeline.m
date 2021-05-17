%% Air Puff
try
    Puff = Event.Puff;
catch
    disp('Old format - no Event structure');
end
for sn=1:nSess
    sess = session_types{use_sessions(sn)};
    
    if isscalar(Puff.start), ps = Puff.start;
    else,                    ps = Puff.start(sn);
    end
    Puff.(sess) = [];
    % get puff times:
    if contains(sess, 'AP') || contains(sess,'puff')
        sess
        tmp = Norm.(sess).t; 
        nTrials = p.nTrials; nTimepoints = size(Norm.(sess).t,1)/nTrials;
        tmp = reshape(tmp,[nTimepoints, nTrials, p.nROIs]);
        
        pufftime.(sess) = ps + arrayfun( @(trial) min(tmp(1,trial,:)), 1:p.nTrials)-min(tmp(:));  %ms
        add_puff = false;
        Puff.(sess) = [pufftime.(sess)' pufftime.(sess)'+Puff.dur];
    end
end
Event.Puff = Puff;

%% Locomotion periods
try
    Loco = Event.Loco;
catch
    disp('Old format - no Event structure');
end
for sn = 1:nSess
    sess = session_types{use_sessions(sn)};
    if ~isempty(Speed.(sess))
    Loco.Detected.(sess) = get_loco_period( Speed.(sess), Loco.use_absolute,    Loco.min_period_separation,   Loco.min_period_size, ...
                                                       [ Loco.on_threshold,   Loco.off_threshold ],         Loco.Avg_window);
    else
        Loco.Detected.(sess)=[];
    end
end

%% Whisking periods
try
    Puff = Event.Whisk;
catch
    disp('Old format - no Event structure');
end

for sn = 1:nSess
    sess = session_types{use_sessions(sn)};
    if exist('VidMI','var')
    Whisk.Detected.(sess) = get_loco_period( VidMI.(sess).(Whisk.vidseg)(:,[2 1]), Whisk.use_absolute,    Whisk.min_period_separation,   Whisk.min_period_size, ...
                                                       [ Whisk.on_threshold,   Whisk.off_threshold ],         Whisk.Avg_window);
    end
end

%% Other Events
try
    for jj=1:length(Event.types)
        fld = Event.types{jj};
        if ~( strcmp(fld, 'Puff') || strcmp(fld, 'Loco') || strcmp(fld, 'Whisk') )
        end
    end
catch
end
%%
save( [exp.analysed, '\Beh.mat'], 'pufftime','Loco','Whisk', 'Puff' ,'-v7.3' );
disp('.....................Saved Behaviour Epochs')