nComp_n = size( PCA.all.soma.proj,1);

% allt = []; maxt=0;
% for jj=1:nSess
%     sess = session_types{use_sessions(jj)};
%     totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
%     allt = [allt; Norm.(sess).t + maxt-totalt];
% end
% 
% maxt=0;
% all_events = Puff;
% concat_all_events.AP = [];
% concat_all_events.puff_out = [];
% concat_all_events.puff_paw = [];
% concat_all_events.whisk = [];
% for jj=1:nSess
% sess = session_types{use_sessions(jj)};
% totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
% concat_all_events.whisk = [concat_all_events.whisk; Whisk.Detected.(sess) + maxt-totalt];
% if contains(sess,'AP')
% all_events.(sess) = Puff.(sess) + maxt-totalt;
% concat_all_events.AP = [concat_all_events.AP; Puff.(sess) + maxt-totalt];
% elseif  contains(sess,'puff_out')
% all_events.(sess) = Puff.(sess) + maxt-totalt;
% concat_all_events.puff_out = [concat_all_events.puff_out; Puff.(sess) + maxt-totalt];
% elseif  contains(sess,'puff_paw')
% all_events.(sess) = Puff.(sess) + maxt-totalt;
% concat_all_events.puff_paw = [concat_all_events.puff_paw; Puff.(sess) + maxt-totalt];
% end
% end

% allspeed = []; maxt=0;
% for jj=1:nSess
%     sess = session_types{use_sessions(jj)};
%     totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
%     temp = Speed.(sess)+[maxt-totalt, 0];
%     allspeed = [allspeed; temp];
% end
% 
% allwhisk = []; maxt=0;
% for jj=1:nSess
%     sess = session_types{use_sessions(jj)};
%     totalt = Norm.(sess).t(end-1); maxt = maxt+totalt;
%     temp = VidMI.(sess).(Whisk.vidseg)+[0, maxt-totalt]; 
%     allwhisk = [allwhisk; temp];
% end



puff_poi = {'AP','whisk_onset'};%,'puff_out','puff_paw'};'AP',
extratime = [500 2000];
EventParams.dsRate = 25;
Beh.loco = beh_all.loco;
Beh.whisk = beh_all.whisk(:,[2,1]);
EventParams.smoothBeh = [100 100];
EventParams.baselineWindow = [-300 -100];

concat_all_events.whisk_onset = [concat_all_events.whisk(:,1),concat_all_events.whisk(:,1)+3000] ;
for pp=1:length(puff_poi)
    puff_sess = puff_poi{pp};
    EventResponses.(puff_sess) = get_event_responses2(  concat_all_events.(puff_sess), activity_all.pc.proj', activity_all.t(:,1), extratime , Beh, EventParams);
    EventResponses.(puff_sess).AlignedT_orig = EventResponses.(puff_sess).AlignedT;
    EventResponses.(puff_sess).AlignedT = EventResponses.(puff_sess).AlignedT+100;

end

% save([exp.analysed,'\EventResponses.mat'],'EventResponses','-v7.3')
PCplus = arrayfun(@(jj) ((prctile(activity_all.pc.eigvec(:,jj),10)>-0.2)*2-1), 1:size(activity_all.pc.proj,1))
fig = plot_event_resp([] ,EventResponses,Puff,PCplus);
% savefig([exp.analysed,'\AP_PCResponses.fig'])

%%

puff_poi = {'whisk_onset','loco_onset'};%,'puff_out','puff_paw'};'AP',,
extratime = [500 2000];
EventParams.dsRate = 25;
Beh.loco = beh_all.loco;
Beh.whisk = beh_all.whisk(:,[2,1]);
EventParams.smoothBeh = [100 100];
EventParams.baselineWindow = [-300 -100];

concat_all_events.whisk_onset = [concat_all_events.whisk(:,1),concat_all_events.whisk(:,1)+3000] ;
concat_all_events.loco_onset = [concat_all_events.loco(:,1),concat_all_events.loco(:,1)+3000] ;
for pp=1:length(puff_poi)
    puff_sess = puff_poi{pp};
    SomaResponses.(puff_sess) = get_event_responses2(  concat_all_events.(puff_sess), activity_all.soma_noPC1', activity_all.t(:,1), extratime , Beh, EventParams);
    SomaResponses.(puff_sess).AlignedT_orig = SomaResponses.(puff_sess).AlignedT;
    SomaResponses.(puff_sess).AlignedT = SomaResponses.(puff_sess).AlignedT+100;

end

save([exp.analysed,'\SomaResponses_noPC1.mat'],'SomaResponses','-v7.3')
fig = plot_event_resp([] ,SomaResponses,Puff,ones(1,p.nROIs));