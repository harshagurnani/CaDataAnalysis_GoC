s_type = {'Moving average', 'Moving median'};
perc = [10, 20, 50];
b_per_trial = [true, false];
c.smooth_scale=20;%ms
for jj=1:2
    for kk = 1:3
        for mm = 1:2
            
            c.smooth_type = s_type{jj};
            c.base_percentile = perc(kk);
            c.baseline_per_trial = b_per_trial(mm);
            tic
            n_new = normalize_and_smooth(r,t,c);
            toc
            figure()
            plot_traces(n_new(:,:,1:30), t(:,:, 1:30), [ 1 3 2])
            title(sprintf('Baseline - \t \t Method: %s, Perc: %d, \t\t Per trial? %d \n Smoothing - \t \t Method: %s, Scale (ms): %d', c.norm_method, c.base_percentile, c.baseline_per_trial, c.smooth_type, c.smooth_scale))
        end
    end
end