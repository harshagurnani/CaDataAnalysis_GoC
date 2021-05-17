%%% mi2 = motion index from whiskerCam
%%% mi1 = motion index from faceCam

%% Starting message
d = dialog('Position',[50 200 420 840],'Name','Testing Video Timing and Syncing', 'Resize','on');

txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 80 380 745],...
           'String',sprintf(['The following script tests basic timing and synchronisation of videos and speed data.\n ---------------- \n', ...
           'OPTION 1: \nEnter the full path+file name (with extension .tif or .txt etc) for all the things you want to test. Also enter the path+name of the file with the timestamps corresponding to the video. Eg: \n',...
           'C:/Harsha/Video/whiskPad.tif  \t OR \nC:/Harsha/Video/WhiskersCam relativetimes.txt.\n',...
           'You can have same body parts or wheel from different videos in the corresponding fields - Different "wheel tifs" in fields Wheel 1 , Wheel 2, Wheel 3 for example.\n\n',...
           'OPTION 2: \nAlternatively, you can give the "all-videos path" and just the filename+extension.\n\n',...
           'OPTION 3: \nFinally, you can give the all-videos path, and just type e/b/w for standard format:.\n',...
           '[Eye/Body/Whiskers]Cam_*.tif \t OR [Eye/Body/Whiskers]Cam-relative times.txt \n',...
           'where * will be wheel/whiskpad/extra1/extra2 depending on category, and the relative times textfile corresponding to the camera will be taken.\n ---------------- \n',...
           'If you have used Image J Macro for stacks, add     0  ( i.e. zero)    instead of actual number of stacks if using e/b/w format \n or      ;0   (i.e. semicolon zero)     after the complete filename if using a non-standard format for filename.\n', ...
           'This 0 imports all files with prefix in the respective field and -stackNN suffix added to the given filename (NN is the substack number). \nDefault value of 600 frames/stack chosen. \n --------------- \n'...
           'If all-videos path is non-empty, it will be APPENDED to all video files/paths (i.e. EXCEPT THE SPEED PATH/FILE)!!!']));

btn = uicontrol('Parent',d,...
           'Position',[135 20 70 25],...
           'String','OK',...
           'Callback','delete(gcf)');

%% Input file paths for analysis      
prompt = {'Wheel 1','Wheel 1 Times','Wheel 2','Wheel 2 Times','Wheel 3','Wheel 3 Times',...
       'Whisker Pad 1','Whisker Pad 1 Times','Whisker Pad 2','Whisker Pad 2 Times',...
       'Extra Area 1.1','Extra Area 1.1 Times', 'Extra Area 1.2', 'Extra Area 1.2 Times',...
       'Extra Area 2.1','Extra Area 2.1 Times', 'Extra Area 2.2','Extra Area 2.2 Times','Speed folder path or path+file', 'All videos folder'};
dlg_title = 'Give full path\filename';
num_lines = 1;
resize='on';

%Default params
defaultans = {['C:\Harsha\wheel.tif'], ['C:\Harsha\BodyCam-relative times.txt'], '','','', '','', '','', '','','','','','','','','','',''};
ans_ROIs = inputdlg(prompt,dlg_title,num_lines,defaultans, resize);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wheel
whl=[];
if ~isempty(ans_ROIs{1}) && ~isempty(ans_ROIs{2})
    whl=[whl,1];
end
if ~isempty(ans_ROIs{3}) && ~isempty(ans_ROIs{4})
    whl=[whl,2];
end
if ~isempty(ans_ROIs{5}) && ~isempty(ans_ROIs{6})
    whl=[whl,3];    
end
wheel = { ans_ROIs{2*whl-1} };
wheel_times = { ans_ROIs{2*whl} };
nWHL = size(whl,2);

% Whisker Pad
whp=[];
if ~isempty(ans_ROIs{7}) && ~isempty(ans_ROIs{8})
    whp=[whp,4];
end
if ~isempty(ans_ROIs{9}) && ~isempty(ans_ROIs{10})
    whp=[whp,5];
end
whiskpad = { ans_ROIs{2*whp-1} };
whiskpad_times = { ans_ROIs{2*whp} };
nWHP = size(whp,2);


% Extra Area 1
ea1=[];
if ~isempty(ans_ROIs{11}) && ~isempty(ans_ROIs{12})
    ea1=[ea1,6];
end
if ~isempty(ans_ROIs{13}) && ~isempty(ans_ROIs{14})
    ea1=[ea1,7];
end
extra1 = { ans_ROIs{2*ea1-1} };
extra1_times = { ans_ROIs{2*ea1} };
nEA1 = size(ea1,2);


% Extra Area 2
ea2=[];
if ~isempty(ans_ROIs{15}) && ~isempty(ans_ROIs{16})
    ea2=[ea2,8];
end
if ~isempty(ans_ROIs{17}) && ~isempty(ans_ROIs{18})
    ea2=[ea2,9];
end
extra2 = { ans_ROIs{2*ea2-1} };
extra2_times = { ans_ROIs{2*ea2} };
nEA2 = size(ea2,2);

% Speed Data
if ~isempty(ans_ROIs{19})
    spd = get_speed_data(ans_ROIs{19});
else
    spd=[];
    speed=[];
end


vd_path = ans_ROIs{20};
if ~isempty(vd_path)
    if ~(vd_path(end) == '\' || vd_path(end) == '/')
        vd_path = [vd_path,'/'];
    end
end
%parpool

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                  Analyse speed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(spd)
    n_trials = size(spd,2);
    speed=[];
    for jj=1:n_trials
        speed=[speed; spd{jj}];

    end
    Last_timestamps.speed = speed(end,1);
    disp(sprintf('\n Loaded speed data. Proceeding to motion index on wheel data... \n'))
else
    disp(sprintf('\n No speed data available. Proceeding to wheel... \n'))
end

% Check last speed timestamp

% Last_timestamps.wheel = cell(size(whl));
% duplications.wheel.timestamps = nan(nWHL,2);
% duplications.wheel.frames = cell(size(whl));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                   Analyse Wheel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nWHL>=1

    mi_wheel = cell(1, nWHL);
    plot_whmi = {};
    
    for jj=1:nWHL
       % Parse standard format for wheel data
       if size(wheel{jj},2) <= 4
           if size(wheel{jj},2) > 1
            nstck_wh = str2double(wheel{jj}(2:end));  
           end
           switch wheel{jj}(1)
               case 'e'
                   wheel{jj} = 'EyeCam_wheel.tif';
               case 'w'
                   wheel{jj} = 'WhiskersCam_wheel.tif';
               case 'b'
                  wheel{jj} = 'BodyCam_wheel.tif';
           end
           
       elseif ~isempty(find(wheel{jj}==';', 1))
           semic = find(wheel{jj}==';', 1);
           time_semic = find(wheel_times{jj}==';',1);
           if ~isempty(time_semic)
               wheel_times{jj} = wheel_times{jj}(1:time_semic-1);
           end
           nstck_wh = str2double(wheel{jj}(semic+1:end));
           wheel{jj} = wheel{jj}(1:semic-1);

       end
       
       %Parse standard format for time data
       if size(wheel_times{jj},2) <= 4
           if size(wheel_times{jj},2) > 1 && isempty(nstck_wh)
            nstck_wh = str2double(wheel_times{jj}(2:end));  
           end
           switch wheel_times{jj}(1)
               case 'e'
                   wheel_times{jj} = 'EyeCam-relative times.txt';
               case 'w'
                   wheel_times{jj} = 'WhiskersCam-relative times.txt';
               case 'b'
                   wheel_times{jj} = 'BodyCam-relative times.txt';
           end
       end
       
       wheel_datatimes{jj} = importdata([vd_path,wheel_times{jj}]);
       Last_timestamps.wheel{jj} = wheel_datatimes{jj}(end,2);
       
       % Number of Duplications in timestamps
       duplications.wheel.timestamps(jj, :) = [size(wheel_datatimes{jj}(:,1),1)-size(unique(wheel_datatimes{jj}(:,1)),1), size(wheel_datatimes{jj}(:,2),1)-size(unique(wheel_datatimes{jj}(:,2)),1)];
       
       disp(sprintf('Computing motion index for wheel %d ...',jj)) %#ok<*DSPS>
       try 
           if ~isempty(nstck_wh)
                    mi_wheel{jj} = simple_motion_index( [vd_path,wheel{jj}], wheel_datatimes{jj} , nstck_wh);
           else 
               mi_wheel{jj} = simple_motion_index( [vd_path,wheel{jj}], wheel_datatimes{jj} );
           end
           
           % Possible duplicated frames: Indices where MI for wheel video jj is zero. 
           duplications.wheel.frames{jj} = find(mi_wheel{jj}(:,1) == 0 );   
       catch ME
            disp(sprintf('Wheel %d MI: %s',jj,ME.message))
       end
       clear nstck_wh
    end

    if ~isempty(speed)
        disp(sprintf('\n Computing correlations between wheel motion as well as speed ...'))
    end
else
    disp(sprintf('\n No wheel data available. Cannot check against speed. \n \n Proceeding to whisker pad data...'))
end


%% Correlations for wheel and speed data
% ----------------------------------------- 

ds_rate = 25;   %Hz - can be higher, if cameras recorded at higher frame rate.
lagbins = 50;   %No. of lag bins on each side of 0


plot_ws={};     % Indices of successful speed-wheel corr
plot_ww={};     % Indices of wheel-wheel corr

%if nWHL>0
for jj=1:size(wheel,2)
    
    if ~isempty(speed)
        try
        [speed_corr{jj}, speed_lag{jj}] = get_lagcorr_ds_synctest( mi_wheel{jj}(:,1), mi_wheel{jj}(:,2),...
                                                            abs(speed(:,2)), speed(:,1), ds_rate, lagbins ) ;  
        %plotws(jj) = 1;%cat(2, plot_ws, num2str(jj));                                                
        plot_ws = cat(2, plot_ws, num2str(jj));
        catch ME
            disp(sprintf('Error with correlation between speed and wheel %d motion index:\n %s', jj, ME.message))
        end
    end
    %temp = zeros(1,nWHL);
    for kk=jj+1:size(wheel,2)
        try
       [ wheel_corr{jj}{kk}, wheel_lag{jj}] = get_lagcorr_ds_synctest( mi_wheel{jj}(:,1), mi_wheel{jj}(:,2),...
                                                                mi_wheel{kk}(:,1), mi_wheel{kk}(:,2),...
                                                                ds_rate, lagbins); 
        %temp(kk) = 1; %cat(2, plot_ww, sprintf('Wheel %d - Wheel %d', jj,kk) );
        plot_ww = cat(2, plot_ww, sprintf('Wheel %d - Wheel %d', jj,kk) );
        catch ME
            disp(sprintf('Error with correlation between wheel %d and wheel %d motion index:\n %s', jj,kk, ME.message))
        end
    end
    %plotww(jj,:)=temp;
end

% for jj=1:size(wheel,2)
%     if plotws(jj)==1
%     plot_ws = cat(2, plot_ws, num2str(jj));  
%     end
%     for kk=jj+1:size(wheel,2)
%         if plotww(jj,kk)==1
%             plot_ww = cat(2, plot_ww, sprintf('Wheel %d - Wheel %d', jj,kk) ); 
%         end
%     end
% end

%end

%% Plot correlation for speed and wheel data
% ----------------------------------------- 

% 1) Lag correlation between speed data and motion index of wheel
if nWHL>=1
    figure;colormap(lines)
    title('Lag correlation between speed data and motion index of wheel')
    for jj=1:nWHL
        try
        plot(speed_lag{jj}, speed_corr{jj})
        catch ME
        end
        hold on
    end
    legend(strcat('Wheel ',plot_ws))
    xlabel( 'Lag (ms)')
    ylabel('Correlation coefficient')
    
    figure;colormap(lines)
    title('Speed data and motion index of wheel')
    plot( speed(:,1), 4*smooth(abs(speed(:,2)),30)/max(abs(speed(:,2)))+0.1 )
    hold on
    for jj=1:nWHL
        try
        plot(mi_wheel{jj}(:,2), mi_wheel{jj}(:,1))
        catch ME
        end
        hold on
    end
    legend(cat(2, 'Speed' ,strcat('Wheel ',plot_ws)))
    xlabel( 'Time (ms)')
    ylabel('Normaised speed or motion index')
end

% 2) Lag correlation between motion indices of different wheel videos
if nWHL>1
    figure;colormap(lines)
    title('Lag correlation between motion indices of different wheel videos')
    for jj=1:nWHL
        for kk=jj+1:nWHL
            try
            plot(wheel_lag{jj}, wheel_corr{jj}{kk})
            catch ME
            end
            hold on
        end
    end
    legend(plot_ww)
    xlabel( 'Lag (ms)')
    ylabel('Correlation coefficient')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%              Analyse Whisker Pad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nWHP>=1
    plot_wpmi={};
    errmsg = cell(1,nWHP);
    mi_whiskpad = cell(1,nWHP);
    for jj=1:nWHP
       if size(whiskpad_times{jj},2) <=4
           if size(whiskpad_times{jj},2) > 1
            nstck_wp = str2double(whiskpad_times{jj}(2:end));  
           end
           switch whiskpad_times{jj}(1)
               case 'e'
                   whiskpad_times{jj} = 'EyeCam-relative times.txt';
                   whiskpad{jj} = 'EyeCam_whiskpad.tif';
               case 'w'
                   whiskpad_times{jj} = 'WhiskersCam-relative times.txt';
                   whiskpad{jj} = 'WhiskersCam_whiskpad.tif';
               case 'b'
                   whiskpad_times{jj} = 'BodyCam-relative times.txt';
                   whiskpad{jj} = 'BodyCam_whiskpad.tif';
           end
        
       elseif ~isempty(find(whiskpad{jj}==';', 1))
           semic = find(whiskpad{jj}==';', 1);
           whiskpad_times{jj} = whiskpad_times{jj}(1:semic-1);
           nstck_wp = whiskpad{jj}(semic+1:end);
           whiskpad{jj} = whiskpad{jj}(1:semic-1);

       end
       whiskpad_datatimes{jj} = importdata([vd_path,whiskpad_times{jj}]);
       Last_timestamps.whiskpad{jj} = whiskpad_datatimes{jj}(end,2);
       % Number of Duplications in timestamps
       duplications.whiskpad.timestamps(jj, 1:2) = [size(whiskpad_datatimes{jj}(:,1),1)-size(unique(whiskpad_datatimes{jj}(:,1)),1), size(whiskpad_datatimes{jj}(:,2),1)-size(unique(whiskpad_datatimes{jj}(:,2)),1)];
       disp(sprintf('Computing motion index for whiskpad %d ...',jj)) %#ok<*DSPS>
       try 
           if ~isempty(nstck_wp)
%                if nstck_wp > 1
                   mi_whiskpad{jj} = simple_motion_index( [vd_path,whiskpad{jj}], whiskpad_datatimes{jj} , nstck_wp);
%                end
           else
               mi_whiskpad{jj} = simple_motion_index( [vd_path,whiskpad{jj}], whiskpad_datatimes{jj});
           end
           % Possible duplicated frames: Indices where MI for whiskpad video jj is zero. 
           duplications.whiskpad.frames{jj} = find(mi_whiskpad{jj}(:,1) == 0 );
           plot_wpmi = cat(2, plot_wpmi, num2str(jj) );
       catch ME
           errmsg{jj}=ME.message;
       end
    end
%     Last_timestamps.whiskpad = ts;
%     duplications.whiskpad.timestamps = ts_dups;
%     duplications.whiskpad.frames = fr_dups;
    for jj=1:nWHP
        if ~isempty(errmsg{jj}) 
            disp(sprintf('Whisker Pad %d MI: %s',jj, errmsg{jj}))
        end
    end

    disp('Computing correlations between whiskpad motion...')
else
    disp(sprintf('\n No whiskpad data available. Proceeding to extra area 1... \n'))
end
% clear ts ts_dups fr_dups errmsg
clear errmsg

%% Correlations for whiskpad motion index
% ----------------------------------------- 

ds_rate = 25;   %Hz - can be higher, if cameras recorded at higher frame rate.
lagbins = 50;   %No. of lag bins on each side of 0


plot_wp={};
%plotwp=zeros(nWHP,nWHP);

% if nWHP>0
for jj=1:nWHP

   %temp = zeros(1,nWHP);
    for kk=jj+1:nWHP
        try
       [ whiskpad_corr{jj}{kk}, whiskpad_lag{jj}] = get_lagcorr_ds_synctest( mi_whiskpad{jj}(:,1), mi_whiskpad{jj}(:,2),...
                                                                mi_whiskpad{kk}(:,1), mi_whiskpad{kk}(:,2),...
                                                                ds_rate, lagbins); 
        %temp(kk) = 1; %cat(2, plot_wp, sprintf('whiskpad %d - whiskpad %d', jj,kk) );
        plot_wp = cat(2, plot_wp, sprintf('Whiskpad %d - Whiskpad %d', jj,kk) ); 
        catch ME
            disp(sprintf('Error with correlation between whiskpad %d and whiskpad %d motion index:\n %s', jj,kk, ME.message))
        end
    end
    %plotwp(jj,:)=temp;
end
% end
% for jj=1:nWHP
% 
%     for kk=jj+1:nWHP
%         if plotwp(jj,kk)==1
%             plot_wp = cat(2, plot_wp, sprintf('Whiskpad %d - Whiskpad %d', jj,kk) ); 
%         end
%     end
% end

%% Plot correlation for whiskpad data
% ----------------------------------------- 

% Lag correlation between motion indices of different whiskpad videos
if nWHP>1
    figure;colormap(lines)
    title('Lag correlation between motion indices of different whiskpad videos')
    for jj=1:nWHP
        for kk=jj+1:nWHP
            try
            plot(whiskpad_lag{jj}, whiskpad_corr{jj}{kk})
            catch ME
            end
            hold on
        end
    end
    legend(plot_wp)
    xlabel( 'Lag (ms)')
    ylabel('Correlation coefficient')
    

end

%Motion indices of different whiskpad videos
if nWHP>=1
    figure;colormap(lines)
    title('Motion indices of different whiskpad videos')
    for jj=1:nWHP

        try
        plot(mi_whiskpad{jj}(:,2), mi_whiskpad{jj}(:,1))
        catch ME
        end
        hold on

    end
    legend(strcat('Whisker Pad ', plot_wpmi))
    xlabel( 'Time (ms)')
    ylabel('Motion index for whisker pad')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Analyse Extra Area 1 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nEA1>=1
%     ts_dups = nan( nEA1, 2);
    errmsg = cell(1,nEA1);
    mi_extra1 = cell(1, nEA1);
    for jj=1:nEA1
        if size(extra1_times{jj},2) <=4
           if size(extra1_times{jj},2) > 1
            nstck = extra1_times{jj}(2);  
           end
           switch extra1_times{jj}(1)
               case 'e'
                   extra1_times{jj} = 'EyeCam-relative times.txt';
                   extra1{jj} = 'EyeCam_extra1.tif';
               case 'w'
                   extra1_times{jj} = 'WhiskersCam-relative times.txt';
                   extra1{jj} = 'WhiskersCam_extra1.tif';
               case 'b'
                   extra1_times{jj} = 'BodyCam-relative times.txt';
                   extra1{jj} = 'BodyCam_extra1.tif';
           end
           
        elseif ~isempty(find(extra1{jj}==';', 1))
           semic = find(extra1{jj}==';', 1);
           extra1_times{jj} = extra1_times{jj}(1:semic-1);
           nstck = extra1{jj}(semic+1:end);
           extra1{jj} = extra1{jj}(1:semic-1);

        end

       extra1_datatimes{jj} = importdata([vd_path,extra1_times{jj}]);
       Last_timestamps.extra1{jj} = extra1_datatimes{jj}(end,2);
       % Number of Duplications in timestamps
       duplications.extra1.timestamps(jj,:) = [size(extra1_datatimes{jj}(:,1),1)-size(unique(extra1_datatimes{jj}(:,1)),1), size(extra1_datatimes{jj}(:,2),1)-size(unique(extra1_datatimes{jj}(:,2)),1)];
       disp(sprintf('Computing motion index for extra area 1 %d ...',jj)) %#ok<*DSPS>
       try 
           if ~isempty(nstck)
%              if nstck>1
                mi_extra1{jj} = simple_motion_index( [vd_path,extra1{jj}], extra1_datatimes{jj}, nstck );
%              end
           else
                mi_extra1{jj} = simple_motion_index( [vd_path,extra1{jj}], extra1_datatimes{jj} );
           end    
           % Possible duplicated frames: Indices where MI for extra1 video jj is zero. 
           duplications.extra1.frames{jj} = find(mi_extra1{jj}(:,1) == 0 );   
       catch ME
           errmsg{jj} = ME.message;
       end
    end
%     Last_timestamps.extra1 = ts;
%     duplications.extra1.timestamps = ts_dups;
%     duplications.extra1.frames = fr_dups;
    for jj=1:nEA1
        if ~isempty(errmsg{jj}) 
            disp(sprintf('Extra 1.%d MI: %s',jj, errmsg{jj}))
        end
    end
    disp('Computing correlations between extra1 motion...')
else
    disp('No extra1 data available. Proceeding to extra area 2...')
end


% clear ts ts_dups fr_dups errmsg
clear errmsg nstck

%% Correlations for extra1 motion index
% ----------------------------------------- 

ds_rate = 25;   %Hz - can be higher, if cameras recorded at higher frame rate.
lagbins = 50;   %No. of lag bins on each side of 0


plot_ea1={};
% plotea1=zeros(nEA1);
for jj=1:nEA1

%    temp = zeros(1, nEA1);
    for kk=jj+1:nEA1
        try
       [ extra1_corr{jj}{kk}, extra1_lag] = get_lagcorr_ds_synctest( mi_extra1{jj}(:,1), mi_extra1{jj}(:,2),...
                                                                mi_extra1{kk}(:,1), mi_extra1{kk}(:,2),...
                                                                ds_rate, lagbins); 
%         temp(kk) = 1;%cat(2, plot_ea1, sprintf('extra1 %d - extra1 %d', jj,kk) );
        plot_ea1 = cat(2, plot_ea1, sprintf('extra1 %d - extra1 %d', jj,kk) );   
        catch ME
            disp(sprintf('Error with correlation between extra1.%d and extra1.%d motion index:\n %s', jj,kk, ME.message))
        end
    end
%     plotea1(jj,:)=temp;
end

% for jj=1:nEA1
% 
%     for kk=jj+1:nEA1
%         if plotea1(jj,kk)==1
%             plot_ea1 = cat(2, plot_ea1, sprintf('extra1 %d - extra1 %d', jj,kk) );      
%         end
%     end
% end


%% Plot correlation for extra1 data
% ----------------------------------------- 

% Lag correlation between motion indices of different extra1 videos
if nEA1>1
    figure;colormap(lines)
    title('Lag correlation between motion indices of different extra1 videos')
    for jj=1:nEA1
        for kk=jj+1:nEA1
            try
            plot(extra1_lag, extra1_corr{jj}{kk})
            catch ME
            end
            hold on
        end
    end
    legend(plot_ea1)
    xlabel( 'Lag (ms)')
    ylabel('Correlation coefficient')
end


% Motion indices of different extra1 videos
if nEA1>=1
    figure;colormap(lines)
    title('Motion indices of different extra1 videos')
    for jj=1:nEA1
 
        try
        plot(mi_extra1{jj}(:,2), mi_extra1{jj}(:,1))
        catch ME
        end
        hold on

    end
    legend(plot_ea1)
    xlabel( 'Time (ms)')
    ylabel('Motion index of extra part 1')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 Analyse Extra Area 2 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nEA2>=1
%     ts_dups = nan(nEA2,2);
    mi_extra2 = cell(1,nEA2);
    errmsg = cell(1,nEA2);
    for jj=1:nEA2
        if size(extra2_times{jj},2) <=4
           if size(extra2_times{jj},2) > 1
            nstck = str2double(extra2_times{jj}(2:end)  );
           end
           switch extra2_times{jj}(1)
               case 'e'
                   extra2_times{jj} = 'EyeCam-relative times.txt';
                   extra2{jj} = 'EyeCam_extra2.tif';
               case 'w'
                   extra2_times{jj} = 'WhiskersCam-relative times.txt';
                   extra2{jj} = 'WhiskersCam_extra2.tif';
               case 'b'
                   extra2_times{jj} = 'BodyCam-relative times.txt';
                   extra2{jj} = 'BodyCam_extra2.tif';
           end
           
       elseif ~isempty(find(extra2{jj}==';', 1))
           semic = find(extra2{jj}==';', 1);
           extra2_times{jj} = extra2_times{jj}(1:semic-1);
           nstck = extra2{jj}(semic+1:end);
           extra2{jj} = extra2{jj}(1:semic-1);

       end
       extra2_datatimes{jj} = importdata([vd_path,extra2_times{jj}]);
       Last_timestamps.extra2{jj} = extra2_datatimes{jj}(end,2);
       % Number of Duplications in timestamps
       duplications.extra2.timestamps(jj, 1:2) = [size(extra2_datatimes{jj}(:,1),1)-size(unique(extra2_datatimes{jj}(:,1)),1), size(extra2_datatimes{jj}(:,2),1)-size(unique(extra2_datatimes{jj}(:,2)),1)];
       disp(sprintf('Computing motion index for extra area 1 %d ...',jj)) %#ok<*DSPS>
       try 
           if ~isempty(nstck)
            mi_extra2{jj} = simple_motion_index( [vd_path,extra2{jj}], extra2_datatimes{jj}, nstck );
           else
            mi_extra2{jj} = simple_motion_index( [vd_path,extra2{jj}], extra2_datatimes{jj} );
           end
           % Possible duplicated frames: Indices where MI for extra2 video jj is zero. 
           duplications.extra2.frames{jj} = find(mi_extra2{jj}(:,1) == 0 );   
       catch ME
           errmsg{jj}=ME.message;
       end
    end
%     Last_timestamps.extra2 = ts;
%     duplications.extra2.timestamps = ts_dups;
%     duplications.extra2.frames = fr_dups;
    for jj=1:nEA2
        if ~isempty(errmsg{jj}) 
            disp(sprintf('Extra 2.%d MI: %s',jj, errmsg{jj}))
        end
    end
    disp('Computing correlations between extra2 motion...')
else
    disp('No extra2 data available.')
end

% clear ts ts_dups fr_dups errmsg
clear errmsg

%% Correlations for extra2 motion index
% ----------------------------------------- 

ds_rate = 25;   %Hz - can be higher, if cameras recorded at higher frame rate.
lagbins = 50;   %No. of lag bins on each side of 0


plot_ea2={};
% plotea2 = zeros(nEA2);
for jj=1:nEA2

%    temp = zeros(1, nEA2);
    for kk=jj+1:nEA2
        try
       [ extra2_corr{jj}{kk}, extra2_lag] = get_lagcorr_ds_synctest( mi_extra2{jj}(:,1), mi_extra2{jj}(:,2),...
                                                                mi_extra2{kk}(:,1), mi_extra2{kk}(:,2),...
                                                                ds_rate, lagbins); 
%         temp = 1; %cat(2, plot_ea2, sprintf('extra2 %d - extra2 %d', jj,kk) );
        plot_ea2 = cat(2, plot_ea2, sprintf('extra2 %d - extra2 %d', jj,kk) ); 
        catch ME
            disp(sprintf('Error with correlation between extra2.%d and extra2.%d motion index:\n %s', jj,kk, ME.message))
        end
    end
%     plotea2(jj,:)=temp;
end
% for jj=1:nEA2
% 
%     for kk=jj+1:nEA2
%         if plotea2(jj,kk)==1
%             plot_ea2 = cat(2, plot_ea2, sprintf('extra2 %d - extra2 %d', jj,kk) ); 
%         end
%     end
% end


%% Plot correlation for extra2 data
% ----------------------------------------- 

% Lag correlation between motion indices of different extra2 videos
if nEA2>1
    figure;colormap(lines)
    title('Lag correlation between motion indices of different extra2 videos')
    for jj=1:nEA2
        for kk=jj+1:nEA2
            try
            plot(extra2_lag, extra2_corr{jj}{kk})
            catch ME
            end
            hold on
        end
    end
    legend(plot_ea2)
    xlabel( 'Lag (ms)')
    ylabel('Correlation coefficient')
end

% Motion indices of different extra2 videos
if nEA2>=1
    figure;colormap(lines)
    title('Motion indices of different extra2 videos')
    for jj=1:nEA2

        try
        plot(mi_extra2{jj}(:,2), mi_extra2{jj}(:,1))
        catch ME
        end
        hold on

    end
    legend(plot_ea2)
    xlabel( 'Time (ms)')
    ylabel('Motion index of extra part 2')
end



%% Analysing video segments



% 
% 
% 
% 
% 
% 
% 
% %%
% last = min( mi1(end), mi2(end) );
% 
% on_times = 0:1000:last-3000;
% off_times = on_times + 3000;
% 
% test_per = [on_times',off_times'];
% 
% un1=
% [crop1, time1] = get_TSOI(mi1(:,1), mi1(:,2), test_per, 0);
% [crop2, time2] = get_TSOI(mi2(:,1), mi2(:,2), test_per, 0);
% 
% 
% for jj=1:60%max(size(crop1))
%     jj
%    [vid_cc{jj},vid_lags{jj} ] = get_lagcorr_ds(crop1{jj}, time1{jj}, crop2{jj}, time2{jj}, 25); 
% end