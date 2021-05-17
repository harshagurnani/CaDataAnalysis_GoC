for id = 1:length(gid)
   roi = gid(id)
   
   
   rawF=[];
   for jj=1:length(use_sessions), sess = session_types{use_sessions(jj)};
   rawF = [rawF; Seg.(sess).r(:,roi)]; end
   rawF = inpaint_nans(rawF,5);
   rawF=rawF/mean(rawF);

   [tau, amp,sigma, events]  = spk_autocalibration( rawF, dt0);
   par.dt = dt0;
   par.tau = tau;
   par.a=amp*1.1;
   par.finetune.sigma = sigma*1.2;
   [spikest fit drift] = spk_est(rawF,par);
   spk_manual(id) = struct('par', par, 'spikes', spikest, 'F', rawF, 'roi', roi, 'fit', fit, 'drift', drift);
   
   hold off
   plot(dt0*[1:length(rawF)]*1e3,smooth(rawF,3)+0.1,'b');
   hold on
   stem(spikest*1e3, -0.15*ones(size(spikest)),'b','Marker','none', 'linewidth',1)
   pause(0.1);
end