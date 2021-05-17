%% Load imaging data
disp('Loading imaging data...')
[p,r,t,n]=load_POI();

%% Load speed data and get locomotion periods
disp('Identifying locomotion periods...')
spd = get_speed_data([p.exp_path, '\Speed_Data\']);

lp=get_loco_period(spd, true, 50, 1000, [0.5 0.02], 30);  %Min period=50 ms, Thresholds, Smoothing window = 30 ms

%% Load whisker and eye cam crop of WhiskerPad, and get motion periods from FaceCam
vid_path=uigetdir(pwd, 'Choose video folder...');
vid_path=[vid_path,'\'];
%vid_path='C:\Users\Harsha\Desktop\videos\14_56_01\';

%WhiskerCam
%%% Uncomment this if whisker cam data is fine.
% disp('Motion index from Whisker Cam...')
% whisk_time= importdata([vid_path, 'WhiskersCam-relative times.txt']);
% mi2 = simple_motion_index([vid_path, 'WhiskersPadCrop.tif'], whisk_time);
% mi2_comp2= motion_map_index2([vid_path, 'WhiskersPadCrop.tif'], whisk_time,0.3,0.9,40);

%EyeCam
disp('Motion index from Whisker pad in Eye Cam...')
eye_time= importdata([vid_path,'EyeCam-relative times.txt']);
mi_eye= simple_motion_index([vid_path,'EyeCam_WPad.tif'], eye_time);
mi_eye_comp2= motion_map_index2([vid_path,'EyeCam_WPad.tif'], eye_time, 0.3,0.9,40);

mi_p =get_loco_period( mi_eye(:,[2 1]), false, 150, 500, [0.1,0.05], 30);


%% Plotting 1

disp('Plotting with locomotion periods...')
shift = 0;
poi_to_plot = min(40, size(n,3));
%%%%%%%%%%%%%%
%%% Figure 1: Patches of locomotion period
figure;
cmap=colormap(lines);
for jj=1:size(lp,1)
hold on
a=area(lp(jj,1:2),repmat(30+poi_to_plot,[1,2]));
a.FaceAlpha=0.25;
a.EdgeAlpha=0;
end

%Plot poi_to_plot ROI
plot_traces(n(:,:,1:poi_to_plot),t(:,:,1:poi_to_plot),[1 3 2])

%Plot motion index
% plot(whisk_time(2:end,2)-shift,smooth(mi2(:,1),20)*4+23.5,'r')
 plot(eye_time(2:end,2)-shift,smooth(mi_eye(:,1),20)*4+5.5+poi_to_plot,'k')

%Plot speed data
for jj=1:p.nTrials
plot(spd{jj}(:,1),smooth(abs(spd{jj}(:,2)),50)+8+poi_to_plot,'Color', cmap(jj,:))
end
xlim([0 spd{p.nTrials}(end,1)])

%% Plotting 2

%%%%%%%%%%%%%%
%%% Figure 2: Patches of Motion index period
figure;
cmap=colormap(lines);
for jj=1:size(mi_p,1)
hold on
a=area(mi_p(jj,1:2)-shift,repmat(30+poi_to_plot,[1,2]));
a.FaceAlpha=0.25;
a.EdgeAlpha=0;
end

%Plot poi_to_plot ROI
plot_traces(n(:,:,1:poi_to_plot),t(:,:,1:poi_to_plot),[1 3 2])

%Plot motion index
%plot(whisk_time(2:end,2)-shift,smooth(mi2(:,1),20)*4+23.5,'r')
plot(eye_time(2:end,2)-shift,smooth(mi_eye(:,1),20)*4+5.5+poi_to_plot,'k')

%Plot speed data
for jj=1:p.nTrials
plot(spd{jj}(:,1),smooth(abs(spd{jj}(:,2)),20)+8+poi_to_plot,'Color', cmap(jj,:))
end
xlim([0 spd{p.nTrials}(end,1)])


%% Correlation analysis
% Breakpoint 

corr_samplerate = 25;   %Hz
corr_maxhalflag = 1000; %ms
corr_padding = [800, 300];   %ms [initial extra, end extra]
disp('Computing lag correlations')
[lags, corrs] = all_period_lagcorr(mi_p, mi_eye(:,1), mi_eye(:,2)-shift,  ... Whisking (MI) periods
                                reshape(n, p.nTimepoints*p.nTrials, p.nROIs), reshape(t,p.nTimepoints*p.nTrials, p.nROIs), ... Fluorescence data
                                corr_samplerate, corr_maxhalflag, corr_padding);

% Combine across ROIs and periods                            
all_rois = permute(squeeze(sum(corrs,1)),[2 1]);    %combine all ROIs for each period
all_lp = squeeze(sum(corrs,3));                     %combine all periods for each ROI

nl=size(lags,2);
ticks=8;
jump=floor(nl/ticks);

figure;
imagesc(all_rois)
set(gca, 'Xtick', 1:jump:nl, 'XTickLabel',lags(1:jump:nl) )
xlabel('Lag (ms)')
ylabel('Whisking period')
title('All ROIs combined')

figure;
imagesc(all_lp)
set(gca, 'Xtick', 1:jump:nl, 'XTickLabel',lags(1:jump:nl))
xlabel('Lag (ms)')
ylabel('ROI')
title('All Periods combined')

figure;plot(lags, sum(all_lp))
title('Summed correlation profile')
xlabel('Lag (ms)')
ylabel('Summed correlation coefficient')