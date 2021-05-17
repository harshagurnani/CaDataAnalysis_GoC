%% Basic data

allData.neurons.f = activity_all.soma;
allData.ROI.keep_f = 1:p.nROIs
allData.neurons.time = activity_all.t;
gid = arrayfun(@(jj) spk_manual(jj).roi, 1:length(spk_manual));
allData.ROI.keep_spikes = gid
allData.neurons.spikes = activity_all.spikes;
allData.neurons.fr = activity_all.FR;
hist(allData.neurons.spikes(:))
allData.params.fr_smooth = 5;

%% 

try
    FOV=p.FOV
catch
    FOV = CorrAnalysis.px2um * 512
end
allData.params.FOV_um=FOV;

%%
ID = strrep([exp.Animal,'__', exp.ID],' ', '')
allData.ID = ID
allEvents.ID = ID
allAnalysed.ID = ID
allData.ROI.Coordinates.X = FOV/512*ROI.Coordinates.ROIs_X
allData.ROI.Coordinates.Y = FOV/512*ROI.Coordinates.ROIs_Y
allData.ROI.Coordinates.Z = ROI.Coordinates.ROIs_Z
allData.ROI.masks = Seg.masks;
allData.ROI.mean_image  = Seg.mean_image;

%% check combined_corr if no beh

plot(allData.neurons.time,allData.neurons.f+[1:p.nROIs])
hold on
plot(beh_all.loco(:,1), beh_all.loco(:,2)/5+p.nROIs+2)
plot(beh_all.whisk(:,2), beh_all.whisk(:,1)*3+p.nROIs+1)
try, plot(beh_all.pupil(:,1), beh_all.pupil(:,2)*3+p.nROIs+3); catch, end


allData.behaviour.speed = beh_all.loco;
allData.behaviour.whiskerMI = beh_all.whisk(:,[2,1]);
allData.raw.wheel = beh_all.encoder;
allData.raw.whiskerMI = beh_all.unsmoothedwhisk(:,[2,1]);
allData.raw.F = activity_all.rawF;

%%
ef=[];
for jj=1:nSess
sess = session_types{use_sessions(jj)};
try
ef = [ef; MC_error.(sess)];
catch
ef = [ef; Seg.(sess).error];
end
end
ef=find(ef)
allData.params.error_frames = ef;
start_t = 1;last_t = 0;
allData.SessionStart.index=[]; allData.SessionStart.time=[];
allData.SessionEnd.index=[]; allData.SessionEnd.time=[];
for jj=1:nSess
sess = session_types{use_sessions(jj)};
allData.SessionStart.index = [allData.SessionStart.index, start_t];
allData.SessionStart.time = [allData.SessionStart.time, allData.neurons.time(start_t, 1)];
seslen = length(Norm.(sess).n);
last_t = last_t + seslen;
start_t = start_t + seslen;
allData.SessionEnd.index = [allData.SessionEnd.index, last_t];
allData.SessionEnd.time = [allData.SessionEnd.time, allData.neurons.time(last_t, 1)];
end
allEvents.SessionStart = allData.SessionStart
allEvents.SessionEnd = allData.SessionEnd
allEvents.SessionEnd
flds = {'ROI_type', 'ROIs', 'trials', 'dwelltime_per_pixel_ms', 'signal_channel', 'smooth_scale', 'maxROI', 'maxTrials', 'nTimepoints', 'acquisition_rate' };
for jj = 1:length(flds), txt = flds{jj};
allData.params.(txt) = p.(txt);
end
allData.params.session_types = session_types(use_sessions)
allData.params.FrameSize = 512

%%
allEvents.locomotion = concat_all_events.loco;
allEvents.whisk = concat_all_events.whisk;
allEvents.puff_whisker = concat_all_events.AP;
allEvents.puff_out = concat_all_events.puff_out;

%%
allEvents.whisk = get_loco_period( allData.behaviour.whiskerMI, Whisk.use_absolute,    Whisk.min_period_separation,   Whisk.min_period_size, ...
[ Whisk.on_threshold,   Whisk.off_threshold ],         Whisk.Avg_window)
figure; hold on
for jj=1:size(allEvents.whisk,1), a=area([allEvents.whisk(jj,1), allEvents.whisk(jj,2)]/1e3, [p.nROIs+5, p.nROIs+5]); a.FaceAlpha = 0.3; end
plot(allData.behaviour.whiskerMI(:,1)/1e3, allData.behaviour.whiskerMI(:,2)*3+p.nROIs+1)
plot(allData.neurons.time/1e3,zscore(allData.neurons.f)/5+[1:p.nROIs])

%%
allEvents.locomotion = get_loco_period( abs(allData.behaviour.speed), Loco.use_absolute,    Loco.min_period_separation,   Loco.min_period_size, ...
[ Loco.on_threshold,   Loco.off_threshold ],         Loco.Avg_window);
figure; hold on
for jj=1:size(allEvents.locomotion,1), a=area([allEvents.locomotion(jj,1), allEvents.locomotion(jj,2)]/1e3, [p.nROIs+5,p.nROIs+5]); a.FaceAlpha = 0.3; end
plot(allData.behaviour.speed(:,1)/1e3, 0.2*abs(allData.behaviour.speed(:,2))+p.nROIs+1)
plot(allData.neurons.time/1e3,zscore(allData.neurons.f)/5+[1:p.nROIs])


%%
save_allData_allAnalysed
allAnalysed.Dimensionality.dff
allAnalysed.Dimensionality.fr