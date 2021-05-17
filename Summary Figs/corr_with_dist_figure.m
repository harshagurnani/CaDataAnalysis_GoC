% Figure for correlations versus density
radius = 20;%um

%Baseline
sess='baseline';
[SortedDist.(sess), DistID.(sess)] = sort(   all_pw_dist.baseline, 'ascend' );
MeanCorr.(sess) = scatstat1( all_pw_dist.baseline, all_baseline_corr, radius, @nanmean );
StdCorr.(sess) = scatstat1( all_pw_dist.baseline, all_baseline_corr, radius, @nanstd );



%Combined
sess='combined';
[SortedDist.(sess), DistID.(sess)] = sort(   all_pw_dist.combined, 'ascend' );
MeanCorr.(sess) = scatstat1( all_pw_dist.combined, all_combined_corr, radius, @nanmean );
StdCorr.(sess) = scatstat1( all_pw_dist.combined, all_combined_corr, radius, @nanstd );



%No PC1
sess='noPC1';
[SortedDist.(sess), DistID.(sess)] = sort(   all_pw_dist.noPC1, 'ascend' );
MeanCorr.(sess) = scatstat1( all_pw_dist.noPC1, all_noPC1_corr, radius, @nanmean );
StdCorr.(sess) = scatstat1( all_pw_dist.noPC1, all_noPC1_corr, radius, @nanstd );

%% plotting
figure;
sess='baseline';
h=fill( [ SortedDist.(sess); flipud(SortedDist.(sess))  ], [MeanCorr.(sess)(DistID.(sess)) - StdCorr.(sess)((DistID.(sess))); flipud( MeanCorr.(sess)(DistID.(sess)) + StdCorr.(sess)((DistID.(sess))) ) ], [0.9 0.9 0.9]);
hold on;
plot( SortedDist.(sess), MeanCorr.(sess)(DistID.(sess)) , 'k','linewidth', 2);
scatter( all_pw_dist.baseline, all_baseline_corr, 2, 'filled')
suptitle(sess)

figure;
sess='combined';
h=fill( [ SortedDist.(sess); flipud(SortedDist.(sess))  ], [MeanCorr.(sess)(DistID.(sess)) - StdCorr.(sess)((DistID.(sess))); flipud( MeanCorr.(sess)(DistID.(sess)) + StdCorr.(sess)((DistID.(sess))) ) ], [0.9 0.9 0.9]);
hold on;
plot( SortedDist.(sess), MeanCorr.(sess)(DistID.(sess)) , 'k','linewidth', 2);
scatter( all_pw_dist.combined, all_combined_corr, 2, 'filled')
suptitle(sess)

figure;
sess='noPC1';
h=fill( [ SortedDist.(sess); flipud(SortedDist.(sess))  ], [MeanCorr.(sess)(DistID.(sess)) - StdCorr.(sess)((DistID.(sess))); flipud( MeanCorr.(sess)(DistID.(sess)) + StdCorr.(sess)((DistID.(sess))) ) ], [0.9 0.9 0.9]);
hold on;
plot( SortedDist.(sess), MeanCorr.(sess)(DistID.(sess)) , 'k','linewidth', 2);
scatter( all_pw_dist.noPC1, all_noPC1_corr, 2, 'filled')
suptitle(sess)
