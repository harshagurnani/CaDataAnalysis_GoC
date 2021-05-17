% Load Whisker tracking data
% file_baseline = 'D:\video tracking\whiskers\hg03_e2\whisk_11_38DeepCut_resnet50_whisk_HG03_e2_1Feb15shuffle1_500000.csv';
% file_AP= 'D:\video tracking\whiskers\hg03_e2\whisk_11_42DeepCut_resnet50_whisk_HG03_e2_1Feb15shuffle1_500000.csv';

n_whisk = 3;
n_poi_per_whisk = 3;
parts_before_whisk = 0;
xy_ids = [2 3]' + (0:n_whisk*n_poi_per_whisk-1)*3 + 3*parts_before_whisk;
lkhood = 4 + (0:n_whisk*n_poi_per_whisk-1)*3 + 3*parts_before_whisk;

% allDLC.baseline = importdata( file_baseline );
% allDLC.AP       = importdata( file_AP );


sess = 'baseline';

%% Position 
WPos.baseline = allDLC.baseline.data(:, xy_ids(:) );
WPos.AP       = allDLC.AP.data(:, xy_ids(:) );

WPos.(sess) = allDLC.(sess).data(:, xy_ids(:) );

PosLk.baseline = allDLC.baseline.data(:, lkhood );
PosLk.AP = allDLC.AP.data(:, lkhood );
PosLk.(sess) = allDLC.(sess).data(:, lkhood );


PosX.(sess) = allDLC.(sess).data(:, xy_ids(1,:) );
PosY.(sess) = allDLC.(sess).data(:, xy_ids(2,:) );

for jj=1:n_poi_per_whisk*n_whisk,
PosX.(sess)(PosLk.(sess)(:,jj)<.9,jj)=NaN;
PosY.(sess)(PosLk.(sess)(:,jj)<.9,jj)=NaN;
PosX.(sess)(PosLk.(sess)(:,jj)<.9,jj) = nanmean(PosX.(sess)(:,jj));%inpaint_nans(PosX.(sess)(:,jj));
PosY.(sess)(PosLk.(sess)(:,jj)<.9,jj) = nanmean(PosY.(sess)(:,jj));%inpaint_nans(PosX.(sess)(:,jj));
end

sess = 'puff_out';

%% whisker angle
%theta - whisker angle
alltheta.(sess) = nan(length(WPos.(sess)), n_whisk);
%method 1: using end points:
for whisker = 1:n_whisk
    id1 = 1+(whisker-1)*n_poi_per_whisk; id2 =  3+(whisker-1)*n_poi_per_whisk;
    dx = PosX.(sess)(:,id1) - PosX.(sess)(:,id2); %allDLC.(sess).data(:, xy_ids(1,id1) ) - allDLC.(sess).data(:, xy_ids(1,id2) );
    dy = PosY.(sess)(:,id1) - PosY.(sess)(:,id2);%allDLC.(sess).data(:, xy_ids(2,id1) ) - allDLC.(sess).data(:, xy_ids(2,id2) );
    alltheta.(sess)(:, whisker) = atan2(dx, dy) * 180/pi;
    shift = pi/4*(1:7);
    ctr = 0;
    while ctr<8 && (min(alltheta.(sess)(:, whisker))<-160 && max(alltheta.(sess)(:, whisker))>160)
        ctr = ctr+1;
        alltheta.(sess)(:, whisker) = mod(shift(ctr)+atan2(dx, dy),2*pi) * 180/pi;
    end
end

%method 2: fit line to all points - get slope and intercept


%% whisker pad
wpad_pts = [0:n_whisk-1]*3 + 1;
wpad.(sess) = [    nanmean(allDLC.(sess).data(:,xy_ids(1,wpad_pts)),2),...
            nanmean(allDLC.(sess).data(:,xy_ids(2,wpad_pts)),2)];
        
%% create smoothing filter
Fs = 300;   %sampling rate (Hz)

%fast_filter
Wp = 45; Ws = 60; %pass and stopband frequencies (Hz) Wp = 5; Ws = 8; 
[filt_fast1, filt_fast2] = get_lp_filter( Fs, Wp, Ws );


%lowpass_filter (setpoint)
Wp = 4; Ws = 6; %pass and stopband frequencies (Hz) %(2,3)
[filt_lp1, filt_lp2] = get_lp_filter( Fs, Wp, Ws );

%slow_filter (for slow state changes)
Wp = 0.5; Ws = 2; %pass and stopband frequencies (Hz)
[filt_slow1, filt_slow2] = get_lp_filter( Fs, Wp, Ws );


%% filter whisking variables


WPos_fast.(sess) = filtfilt( filt_fast1, filt_fast2, WPos.(sess) );
WPos_lp.(sess) = filtfilt( filt_lp1, filt_lp2, WPos.(sess) );
WPos_slow.(sess) = filtfilt( filt_slow1, filt_slow2, WPos.(sess) );


alltheta_fast.(sess) = filtfilt( filt_fast1, filt_fast2, alltheta.(sess) );
alltheta_lp.(sess) = filtfilt( filt_lp1, filt_lp2, alltheta.(sess) );
alltheta_slow.(sess) = filtfilt( filt_slow1, filt_slow2, alltheta.(sess) );


wpad_fast.(sess) = filtfilt( filt_fast1, filt_fast2, wpad.(sess) );
wpad_lp.(sess) = filtfilt( filt_lp1, filt_lp2, wpad.(sess) );
wpad_slow.(sess) = filtfilt( filt_slow1, filt_slow2, wpad.(sess) );

VidMI_fast.(sess) = filtfilt( filt_fast1, filt_fast2, [0;VidMI.(sess).whiskpadR_wcam(:,1)]);
VidMI_lp.(sess) = filtfilt( filt_lp1, filt_lp2, [0;VidMI.(sess).whiskpadR_wcam(:,1)] );
VidMI_slow.(sess) = filtfilt( filt_slow1, filt_slow2, [0;VidMI.(sess).whiskpadR_wcam(:,1)] );


WhiskData.(sess) = struct('WPos_fast', WPos_fast.(sess), 'WPos_lp', WPos_lp.(sess), 'WPos_slow', WPos_slow.(sess), ...
                          'alltheta_fast', alltheta_fast.(sess), 'alltheta_lp', alltheta_lp.(sess), 'alltheta_slow', alltheta_slow.(sess), ...
                          'wpad_fast', wpad_fast.(sess), 'wpad_lp', wpad_lp.(sess), 'wpad_slow', wpad_slow.(sess), ...
                          'VidMI_fast', VidMI_fast.(sess), 'VidMI_lp', VidMI_lp.(sess), 'VidMI_slow', VidMI_slow.(sess),...
                          'VidMI', [0; VidMI.(sess).whiskpadR_wcam(:,1)]);
                      
%% Amplitude of whisking
WhiskAmp.(sess) = nan(size(alltheta_lp.(sess)));
WhiskHil.(sess) = WhiskAmp.(sess);
for whisk = 1:n_whisk
	WhiskHil.(sess)(:,whisk) = hilbert(alltheta_fast.(sess)(:,whisk)-alltheta_lp.(sess)(:,whisk));
    WhiskAmp.(sess)(:,whisk) = abs( WhiskHil.(sess)(:,whisk) );
end