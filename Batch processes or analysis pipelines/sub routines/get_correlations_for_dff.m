%% Pairwise Correlations
%Total correlations
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   all_n = Norm.(sess).n;
   if any(isnan(Norm.(sess).n)) % remove NaN
       for jj=1:p.nROIs, all_n(:,jj) = inpaint_nans(Norm.(sess).n(:,jj)); end
   end
   
   [Corr.pw.(sess).all.lags, Corr.pw.(sess).all.corr, ~, Corr.pw.(sess).all.conf_interval, Corr.pw.(sess).all.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( [0 totalT], all_n, Norm.(sess).t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
end
% 
% % combined traces!
Norm.all.n = [];
Norm.all.t = [];
maxT = 0; 
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   maxT = maxT + totalT;
   all_n = Norm.(sess).n;
   if any(isnan(Norm.(sess).n)) % remove NaN
       for jj=1:p.nROIs, all_n(:,jj) = inpaint_nans(Norm.(sess).n(:,jj)); end
   end
   Norm.all.n = [Norm.all.n; all_n];
   Norm.all.t = [Norm.all.t; Norm.(sess).t+ maxT-totalT+(sn-1)*Corr.pw.params.interSess];
end
[Corr.pw.combined.all.lags, Corr.pw.combined.all.corr, ~, Corr.pw.combined.all.conf_interval, Corr.pw.combined.all.corr_sig ] = ...
                        all_period_lagcorr_pairwise_dff_new( [0 maxT], Norm.all.n, Norm.all.t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );

%% during beh                 
% During running
if exist('Loco', 'var')
    if isfield(Loco, 'Detected')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   if~isempty(Loco.Detected.(sess))
   [Corr.pw.(sess).loco.lags, Corr.pw.(sess).loco.corr, Corr.pw.(sess).loco.period, Corr.pw.(sess).loco.conf_interval, Corr.pw.(sess).loco.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( Loco.Detected.(sess), Norm.(sess).n, Norm.(sess).t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.loco_maxlag, Corr.pw.params.loco.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.loco_doconcat );
   end
end
end
end

% During whisking
if exist('Whisk', 'var')
if isfield(Whisk, 'Detected')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   [Corr.pw.(sess).whisk.lags, Corr.pw.(sess).whisk.corr, Corr.pw.(sess).whisk.period, Corr.pw.(sess).whisk.conf_interval, Corr.pw.(sess).whisk.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( Whisk.Detected.(sess), Norm.(sess).n, Norm.(sess).t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.whisk_maxlag, Corr.pw.params.whisk.extratime, [], Corr.pw.params.nshuffle, false, Corr.pw.params.whisk_doconcat );
end
end
end


% During air puff
if exist( 'Puff', 'var')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   if contains(sess, 'AP') || contains(sess,'puff' ) && isfield( Puff, sess)
   [Corr.pw.(sess).puff.lags, Corr.pw.(sess).puff.corr, Corr.pw.(sess).puff.period, Corr.pw.(sess).puff.conf_interval, Corr.pw.(sess).puff.corr_sig ] =  ...
                        all_period_lagcorr_pairwise_dff_new( Puff.(sess), Norm.(sess).n, Norm.(sess).t, ...
                                Corr.pw.params.ds_rate, Corr.pw.params.puff_maxlag, Corr.pw.params.puff.extratime, [], Corr.pw.params.nshuffle, false );
   end
end
end

% % Background cross correlations
% for sn = 1:nSess
%    sess = session_types{use_sessions(sn)};
%    totalT = Norm.(sess).t(end);
%    [Corr.pw.(sess).background.lags, Corr.pw.(sess).background.corr, ~, Corr.pw.(sess).background.conf_interval, Corr.pw.(sess).background.corr_sig ] = ...
%                         all_period_lagcorr_pairwise_dff_new( [0 totalT], Norm.(sess).background, Norm.(sess).t, ...
%                                 Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false, Corr.pw.params.puff_doconcat );
% end
%---------------------------------------------------

%% Behavioural Correlations
%--------------------------------------------------
% With running
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   
   if ~isempty(Speed.(sess))
   
   loco         = smooth(-Speed.(sess)(:,2),50);                           % inverted encoder trace
   runspeed     = -Speed.(sess)(:,2); runspeed(runspeed<0) = 0;            % only forward run speed
   runspeed     = smooth(runspeed, 50);
   movement     = smooth(abs(Speed.(sess)(:,2)),50);                       % any movement rectified
   
   
   [Corr.beh.(sess).loco.lags, Corr.beh.(sess).loco.corr, Corr.beh.(sess).loco.period, Corr.beh.(sess).loco.conf_interval, Corr.beh.(sess).loco.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], loco, Speed.(sess)(:,1), ...
                                                            Norm.(sess).n, Norm.(sess).t, ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
    [Corr.beh.(sess).runspeed.lags, Corr.beh.(sess).runspeed.corr, Corr.beh.(sess).runspeed.period, Corr.beh.(sess).runspeed.conf_interval, Corr.beh.(sess).runspeed.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], runspeed, Speed.(sess)(:,1), ...
                                                            Norm.(sess).n, Norm.(sess).t, ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
    [Corr.beh.(sess).mov.lags, Corr.beh.(sess).mov.corr, Corr.beh.(sess).mov.period, Corr.beh.(sess).mov.conf_interval, Corr.beh.(sess).mov.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], movement, Speed.(sess)(:,1), ...
                                                            Norm.(sess).n, Norm.(sess).t, ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
   end                            
%    [Corr.beh.(sess).background.loco.lags, Corr.beh.(sess).background.loco.corr, Corr.beh.(sess).background.loco.period, Corr.beh.(sess).background.loco.conf_interval, Corr.beh.(sess).background.loco.corr_sig ] =  ...
%                         all_period_lagcorr_beh_and_dff_new( [0 totalT], loco, Speed.(sess)(:,1), ...
%                                                             Norm.(sess).background, Norm.(sess).t, ...
%                                                             Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );                                                     
%     [Corr.beh.(sess).background.runspeed.lags, Corr.beh.(sess).background.runspeed.corr, Corr.beh.(sess).background.runspeed.period, Corr.beh.(sess).background.runspeed.conf_interval, Corr.beh.(sess).background.runspeed.corr_sig ] =  ...
%                         all_period_lagcorr_beh_and_dff_new( [0 totalT], runspeed, Speed.(sess)(:,1), ...
%                                                             Norm.(sess).background, Norm.(sess).t, ...
%                                                             Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );                                                     
%     [Corr.beh.(sess).background.mov.lags, Corr.beh.(sess).background.mov.corr, Corr.beh.(sess).background.mov.period, Corr.beh.(sess).background.mov.conf_interval, Corr.beh.(sess).background.mov.corr_sig ] =  ...
%                         all_period_lagcorr_beh_and_dff_new( [0 totalT], movement, Speed.(sess)(:,1), ...
%                                                             Norm.(sess).background, Norm.(sess).t, ...
%                                                             Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );                                                     
%                                                         
end

% With Whisking
if exist('VidMI', 'var') && exist('Whisk', 'var')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   tbin = floor(10/nanmean(diff(VidMI.(sess).(Whisk.vidseg)(:,2))));
   [~,ia,~] = unique(VidMI.(sess).(Whisk.vidseg)(:,2));
   uniqueWhisk = smooth(VidMI.(sess).(Whisk.vidseg),tbin);
   uniqueWhisk = uniqueWhisk(ia,1);
   uniqueWhiskTime = VidMI.(sess).(Whisk.vidseg)(ia,2);
   
   [Corr.beh.(sess).whisk.lags, Corr.beh.(sess).whisk.corr, Corr.beh.(sess).whisk.period, Corr.beh.(sess).whisk.conf_interval, Corr.beh.(sess).whisk.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], uniqueWhisk, uniqueWhiskTime,...
                                                            Norm.(sess).n, Norm.(sess).t,  ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
%    [Corr.beh.(sess).background.whisk.lags, Corr.beh.(sess).background.whisk.corr, Corr.beh.(sess).background.whisk.period, Corr.beh.(sess).background.whisk.conf_interval, Corr.beh.(sess).background.whisk.corr_sig ] =  ...
%                         all_period_lagcorr_beh_and_dff_new( [0 totalT], uniqueWhisk, uniqueWhiskTime,...
%                                                             Norm.(sess).background, Norm.(sess).t,  ...
%                                                             Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
end
end

% With Pupil
if exist('Pupil','var')
for sn = 1:nSess
   sess = session_types{use_sessions(sn)};
   totalT = Norm.(sess).t(end);
   [Corr.beh.(sess).pupil.lags, Corr.beh.(sess).pupil.corr, Corr.beh.(sess).pupil.period, Corr.beh.(sess).pupil.conf_interval, Corr.beh.(sess).pupil.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], Pupil.(sess).longAxis,Vidtime.(sess).ecam(:,2),...
                                                            Norm.(sess).n, Norm.(sess).t,  ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
   [Corr.beh.(sess).background.pupil.lags, Corr.beh.(sess).background.pupil.corr, Corr.beh.(sess).background.pupil.period, Corr.beh.(sess).background.pupil.conf_interval, Corr.beh.(sess).background.pupil.corr_sig ] =  ...
                        all_period_lagcorr_beh_and_dff_new( [0 totalT], Pupil.(sess).longAxis,Vidtime.(sess).ecam(:,2),...
                                                            Norm.(sess).background, Norm.(sess).t,  ...
                                                            Corr.beh.params.ds_rate, Corr.beh.params.maxlag, [], [], Corr.beh.params.nshuffle, false );
end    
    
end



save([exp.analysed, '\correlations.mat'], 'Corr','p','-v7.3');
