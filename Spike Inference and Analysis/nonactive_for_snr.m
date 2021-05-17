if ~exist('tm','var')
    % Get no whisk Period
    Whisk.use_absolute          = true;
    Whisk.min_period_separation = 500;    % ms
    Whisk.min_period_size       = 1000;   % ms
    Whisk.on_threshold          = 0.09;
    Whisk.off_threshold         = 0.1;
    Whisk.Avg_window            = 100;    % ms

    wper=get_loco_period( beh_all.whisk(:,[2,1]), Whisk.use_absolute,    Whisk.min_period_separation,   Whisk.min_period_size, ...
    [ Whisk.on_threshold*0.85,   Whisk.off_threshold ],         Whisk.Avg_window);

    nROI = length(gid);

    t = activity_all.t(:,1);

    tm = activity_all.t(:,1);
    for jj=1:length(wper),
    start = max(1,find(t>wper(jj,1)-1500,1)-1);
    last = find(t>wper(jj,2)+1500,1); if isnan(last), last = length(activity_all.t); end
    tm(start:last)=NaN;
    end
    
end

usetm = find(~isnan(tm));

% Sigma to use for redoing spike inference
snr=arrayfun(@(jj) std(spk_manual(jj).F(usetm)),1:nROI);