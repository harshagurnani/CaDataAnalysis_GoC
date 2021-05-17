rawF=[];
for jj=1:length(use_sessions), sess = session_types{use_sessions(jj)};
rawF = [rawF; Seg.(sess).r(:,roi)]; end
rawF = inpaint_nans(rawF,5);
rawF=rawF-min(rawF);
rawF=rawF/max(rawF);
