% get concatenated behaviour epochs
interSess = 1000;
maxt=0;
beh_all.wheel = [];

% for jj=1:nSess
% sess = session_types{use_sessions(jj)};
% totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;
% 
% 
% for jj=1:p.nROIs,
% activity_all.t(:,jj)=inpaint_nans(activity_all.t(:,jj));end
% % concatenate behaviour
% [~,id]=unique(VidMI.(sess).limbwheel_ecam(:,2));
% id = sort(id,'ascend');
% uwhisk=smooth(VidMI.(sess).limbwheel_ecam(:,1),3);
% beh_all.wheel = [beh_all.wheel; [VidMI.(sess).limbwheel_ecam(:,2),uwhisk] + [maxt-totalt, 0] ];
% end
%%
concat_all_events.AP = [];
concat_all_events.puff_out = [];
concat_all_events.puff_paw = [];
concat_all_events.contra_puff = [];
concat_all_events.puff_paw = [];
concat_all_events.whisk = [];
concat_all_events.loco =[];

activity_all.soma = [];
activity_all.t = [];
activity_all.rawF =[];
% activity_all.pc= [];

beh_all.loco = [];
beh_all.whisk = [];
beh_all.unsmoothedwhisk = [];
beh_all.encoder = [];

for jj=1:nSess
sess = session_types{use_sessions(jj)}
totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;

% concatenate behaviour
try,concat_all_events.whisk = [concat_all_events.whisk; Whisk.Detected.(sess) + maxt-totalt]; catch, end
concat_all_events.loco  = [concat_all_events.loco;  Loco.Detected.(sess)  + maxt-totalt];

try,
if contains(sess,'AP')
concat_all_events.AP = [concat_all_events.AP; Puff.(sess) + maxt-totalt];
elseif  contains(sess,'puff_out')
concat_all_events.puff_out = [concat_all_events.puff_out; Puff.(sess) + maxt-totalt];
elseif  contains(sess,'puff_paw')
concat_all_events.puff_paw = [concat_all_events.puff_paw; Puff.(sess) + maxt-totalt];
elseif  contains(sess,'contra_puff')
concat_all_events.contra_puff = [concat_all_events.contra_puff; Puff.(sess) + maxt-totalt];
end
catch,
end
% 
% 
% concatenate activity
activity_all.soma = [activity_all.soma; Norm.(sess).n];
activity_all.rawF = [activity_all.rawF; Seg.(sess).r];
activity_all.t    = [activity_all.t; Norm.(sess).t + maxt-totalt];

for kk=1:p.nROIs, activity_all.soma(:,kk)=inpaint_nans(activity_all.soma(:,kk));
activity_all.t(:,kk)=inpaint_nans(activity_all.t(:,kk));end
% concatenate behaviour
if ~isempty(Speed.(sess))
beh_all.loco = [beh_all.loco; [Speed.(sess)(:,1),smooth(-Speed.(sess)(:,2),50)] + [maxt-totalt, 0] ];
beh_all.encoder = [beh_all.encoder; [Speed.(sess)(:,1),-Speed.(sess)(:,2)] + [maxt-totalt, 0] ];
end

try
[~,id]=unique(VidMI.(sess).(Whisk.vidseg)(:,2));
id = sort(id,'ascend');
uwhisk=smooth(VidMI.(sess).(Whisk.vidseg)(:,1),30);
beh_all.whisk = [beh_all.whisk; [uwhisk(id),  VidMI.(sess).(Whisk.vidseg)(id,2)]+ [0 maxt-totalt] ];
beh_all.unsmoothedwhisk = [beh_all.unsmoothedwhisk; [VidMI.(sess).(Whisk.vidseg)(id,1),  VidMI.(sess).(Whisk.vidseg)(id,2)]+ [0 maxt-totalt] ];
catch, end
% [~,id]=unique(VidMI.(sess).limbwheel_ecam(:,2));
% id = sort(id,'ascend');
% uwhisk=smooth(VidMI.(sess).limbwheel_ecam(:,1),3);
% beh_all.wheel = [beh_all.wheel; [VidMI.(sess).limbwheel_ecam(:,2),uwhisk] + [maxt-totalt, 0] ];
% 
end
activity_all.pc= test_pca_on_fr(activity_all.soma,activity_all.t);
% save([exp.analysed, '\conactenated.mat'], 'beh_all','activity_all','concat_all_events','-v7.3');

%% pairwise correlations 
% get combined correlations
% All
[Corr.pw.combined.all.lags, Corr.pw.combined.all.corr, ~, Corr.pw.combined.all.conf_interval, Corr.pw.combined.all.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxt], activity_all.soma, activity_all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
[PC_Corr.combined.pw.all.lags, PC_Corr.combined.pw.all.corr, ~, PC_Corr.combined.pw.all.conf_interval, PC_Corr.combined.pw.all.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxt], activity_all.pc.proj', repmat(activity_all.t(:,1),1,size(activity_all.pc.proj',2)), ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
% During running
if~isempty(concat_all_events.loco)
   [Corr.pw.combined.loco.lags, Corr.pw.combined.loco.corr, Corr.pw.combined.loco.period, Corr.pw.combined.loco.conf_interval, Corr.pw.combined.loco.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( concat_all_events.loco, activity_all.soma, activity_all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.loco_maxlag, Corr.pw.params.loco.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.loco_doconcat );
   [PC_Corr.combined.pw.loco.lags, PC_Corr.combined.pw.loco.corr, PC_Corr.combined.pw.loco.period, PC_Corr.combined.pw.loco.conf_interval, PC_Corr.combined.pw.loco.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( concat_all_events.loco, activity_all.pc.proj', repmat(activity_all.t(:,1),1,size(activity_all.pc.proj',2)), ...
                                Corr.pw.params.ds_rate, Corr.pw.params.loco_maxlag, Corr.pw.params.loco.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.loco_doconcat );

end

% During whisking
if~isempty(concat_all_events.whisk)
   [Corr.pw.combined.whisk.lags, Corr.pw.combined.whisk.corr, Corr.pw.combined.whisk.period, Corr.pw.combined.whisk.conf_interval, Corr.pw.combined.whisk.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( concat_all_events.whisk, activity_all.soma, activity_all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.whisk_maxlag, Corr.pw.params.whisk.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.whisk_doconcat );
   [PC_Corr.combined.pw.whisk.lags, PC_Corr.combined.pw.whisk.corr, PC_Corr.combined.pw.whisk.period, PC_Corr.combined.pw.whisk.conf_interval, PC_Corr.combined.pw.whisk.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( concat_all_events.whisk, activity_all.pc.proj', repmat(activity_all.t(:,1),1,size(activity_all.pc.proj',2)), ...
                                Corr.pw.params.ds_rate, Corr.pw.params.whisk_maxlag, Corr.pw.params.whisk.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.whisk_doconcat );

end


%% correlations with behaviour
[Corr.combined.beh.loco.lags, Corr.combined.beh.loco.corr, Corr.combined.beh.loco.period, Corr.combined.beh.loco.conf_interval, Corr.combined.beh.loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.loco(:,2), beh_all.loco(:,1), ...
                                                             activity_all.soma, repmat(activity_all.t(:,1),1,p.nROIs), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );

[PC_Corr.combined.beh.loco.lags, PC_Corr.combined.beh.loco.corr, PC_Corr.combined.beh.loco.period, PC_Corr.combined.beh.loco.conf_interval, PC_Corr.combined.beh.loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.loco(:,2), beh_all.loco(:,1), ...
                                                            activity_all.pc.proj', repmat(activity_all.t(:,1),1,size(activity_all.pc.proj',2)), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );

if ~isempty( beh_all.whisk)
[Corr.combined.beh.whisk.lags, Corr.combined.beh.whisk.corr, Corr.combined.beh.whisk.period, Corr.combined.beh.whisk.conf_interval, Corr.combined.beh.whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.whisk(:,1), beh_all.whisk(:,2), ...
                                                             activity_all.soma, repmat(activity_all.t(:,1),1,p.nROIs), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );


[PC_Corr.combined.beh.whisk.lags, PC_Corr.combined.beh.whisk.corr, PC_Corr.combined.beh.whisk.period, PC_Corr.combined.beh.whisk.conf_interval, PC_Corr.combined.beh.whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.whisk(:,1), beh_all.whisk(:,2), ...
                                                            activity_all.pc.proj', repmat(activity_all.t(:,1),1,size(activity_all.pc.proj',2)), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );
end

%% Corr without top component
activity_all.soma_noPC1 = subspace_svd( activity_all.soma', -1 );
figure;
plot(activity_all.t(:,1), activity_all.soma_noPC1' + [1:p.nROIs]*.4); hold on
plot(beh_all.loco(:,1), beh_all.loco(:,2)/15); 
%plot(beh_all.whisk(:,2), beh_all.whisk(:,1)*5-3)

%pw corr
[Corr.pw.noPC1.all.lags, Corr.pw.noPC1.all.corr, ~, Corr.pw.noPC1.all.conf_interval, Corr.pw.noPC1.all.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxt], activity_all.soma_noPC1', activity_all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
if ~isempty( concat_all_events.whisk)
[Corr.pw.noPC1.whisk.lags, Corr.pw.noPC1.whisk.corr, ~, Corr.pw.noPC1.whisk.conf_interval, Corr.pw.noPC1.whisk.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( concat_all_events.whisk, activity_all.soma_noPC1', activity_all.t, ...
                               Corr.pw.params.ds_rate, Corr.pw.params.whisk_maxlag, Corr.pw.params.whisk.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.whisk_doconcat  );
end
%beh corr
[Corr.beh.noPC1.loco.lags, Corr.beh.noPC1.loco.corr, Corr.beh.noPC1.loco.period, Corr.beh.noPC1.loco.conf_interval, Corr.beh.noPC1.loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.loco(:,2), beh_all.loco(:,1), ...
                                                            activity_all.soma_noPC1', repmat(activity_all.t(:,1),1,p.nROIs), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );

if ~isempty(beh_all.whisk)
[Corr.beh.noPC1.whisk.lags, Corr.beh.noPC1.whisk.corr, Corr.beh.noPC1.whisk.period, Corr.beh.noPC1.whisk.conf_interval, Corr.beh.noPC1.whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 maxt], beh_all.whisk(:,1), beh_all.whisk(:,2), ...
                                                            activity_all.soma_noPC1', repmat(activity_all.t(:,1),1,p.nROIs), ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false  );
end
%%
save([exp.analysed, '\conactenated.mat'], 'beh_all','activity_all','concat_all_events','-v7.3');
save([exp.analysed, '\PC_Corr.mat'], 'PC_Corr','-v7.3');
save([exp.analysed, '\correlations.mat'], 'Corr','p','-v7.3');

