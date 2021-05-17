id=id+1
roi=gid(id)
rawF=[];
for jj=1:length(use_sessions), sess = session_types{use_sessions(jj)};
rawF = [rawF; Seg.(sess).r(:,roi)]; end
rawF = inpaint_nans(rawF,5);
rawF=rawF-min(rawF);
rawF=rawF/max(rawF);
hold off;
plot(dt0*[1:length(rawF)]*1e3,smooth(rawF,3)+0.1,'b'); hold on
stem(spk_manual(id).spikes*1e3, -0.15*ones(length(spk_manual(id).spikes),1),'b','Marker','none', 'linewidth',1)


%to rerun = id=
%to remove id=