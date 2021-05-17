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
indx.whl=[];
if ~isempty(ans_ROIs{1}) && ~isempty(ans_ROIs{2})
    indx.whl=[indx.whl,1];
end
if ~isempty(ans_ROIs{3}) && ~isempty(ans_ROIs{4})
    indx.whl=[indx.whl,2];
end
if ~isempty(ans_ROIs{5}) && ~isempty(ans_ROIs{6})
    indx.whl=[indx.whl,3];    
end
tifnames.wheel = { ans_ROIs{2*indx.whl-1} };
timefiles.wheel = { ans_ROIs{2*indx.whl} };
% nindx.whl = size(indx.whl,2);

% Whisker Pad
indx.indx.whp=[];
if ~isempty(ans_ROIs{7}) && ~isempty(ans_ROIs{8})
    indx.indx.whp=[indx.indx.whp,4];
end
if ~isempty(ans_ROIs{9}) && ~isempty(ans_ROIs{10})
    indx.indx.whp=[indx.indx.whp,5];
end
tifnames.whiskpad = { ans_ROIs{2*indx.indx.whp-1} };
timefiles.whiskpad = { ans_ROIs{2*indx.indx.whp} };
% nindx.indx.whp = size(indx.indx.whp,2);


% Extra Area 1
indx.ea1=[];
if ~isempty(ans_ROIs{11}) && ~isempty(ans_ROIs{12})
    indx.ea1=[indx.ea1,6];
end
if ~isempty(ans_ROIs{13}) && ~isempty(ans_ROIs{14})
    indx.ea1=[indx.ea1,7];
end
tifnames.extra1 = { ans_ROIs{2*indx.ea1-1} };
timefiles.extra1 = { ans_ROIs{2*indx.ea1} };
% nEA1 = size(ea1,2);


% Extra Area 2
indx.ea2=[];
if ~isempty(ans_ROIs{15}) && ~isempty(ans_ROIs{16})
    indx.ea2=[indx.ea2,8];
end
if ~isempty(ans_ROIs{17}) && ~isempty(ans_ROIs{18})
    indx.ea2=[indx.ea2,9];
end
tifnames.extra2 = { ans_ROIs{2*indx.ea2-1} };
timefiles.extra2 = { ans_ROIs{2*indx.ea2} };
% nEA2 = size(ea2,2);

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
    disp(sprintf('\n Loaded speed data. Proceeding to motion index on wheel data... \n')) %#ok<*DSPS>
else
    disp(sprintf('\n No speed data available. Proceeding to wheel... \n'))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Analyse Wheel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Last_timestamp.wheel, duplications.wheel, Wheel.MI, Wheel.RefDist, Wheel.CrossCorr, Wheel.CrossCorrLag, Wheel.SpeedCorr, Wheel.SpeedCorrLag ] = process_all_copies( indx.whl, vd_path, tifnames.wheel, timefiles.wheel, 'wheel', speed );
%Wheel.SpeedCorr, Wheel.SpeedCorrLag
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Analyse Whisker Pad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Last_timestamp.whiskpad, duplications.whiskpad, whiskpad.MI, whiskpad.RefDist, whiskpad.CrossCorr, whiskpad.CrossCorrLag,~, ~ ] = process_all_copies( indx.indx.whp, vd_path, tifnames.whiskpad, timefiles.whiskpad, 'whiskpad' );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Analyse Extra Area 1 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Last_timestamp.extra1, duplications.extra1, extra1.MI, extra1.RefDist, extra1.CrossCorr, extra1.CrossCorrLag, ~, ~] = process_all_copies( indx.ea1, vd_path, tifnames.extra1, timefiles.extra1, 'extra1' );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Analyse Extra Area 2 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Last_timestamp.extra2, duplications.extra2, extra2.MI, extra2.RefDist, extra2.CrossCorr, extra2.CrossCorrLag, ~, ~ ] = process_all_copies( indx.ea2, vd_path, tifnames.extra2, timefiles.extra2, 'extra2' );


clear defaultans resize num_lines dlg_title prompt
