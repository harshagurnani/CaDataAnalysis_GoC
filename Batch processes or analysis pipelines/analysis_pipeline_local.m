clear variables
% NOTE - 
% You can adjust different settings, but code within dotted lines should
% not be changed unless necessary as it is mostly automated computations


% << SET >>

%% What Do You Want

% what measures?        
get_dff     = true;                                                         % (MC/Segmentation for patch data +) Normalisation to get dF/F
get_speed   = true;                                                         % Read wheel encoder data
get_videoMI = false;                                                         % Use if video_parts are cropped -> get Motion Index from Videos

get_vidpca  = false;                                                        % Dim reduction of motion index/behaviour using PCA

% what analyses
get_beh_epoch = true;                                                       % Different Behavioural Epochs
get_corr    = true;                                                         % Pairwise and behavioural correlations of dF/F
get_pca     = true;                                                         % PCA of dF/F (and background)

% what plots?
get_plot_simple = true;                                                     % Simple plot - (only mean patch,) dff traces and speed
get_plot_beh    = true;                                                     % can add other beh info available - like mi, vid pca etc
get_puff        = true;                                                     % add air puffs

% what sessions?
session_types   = {'baseline', 'AP', 'AP2', ...                             %these will be used as suffix identifiers
                   'baseline2', 'puff_out', 'puff_paw', 'baseline3', ...
                   'contra_puff', 'extra2', 'extra3' };                          % extras are buffers for additional session types - or different baselines. You can add as many additional session types as you want.
                                                                            % Rename them to relevant descriptor.                                                                           
use_sessions    = [1 2 3];% 2  3 4 5 6 7];                                     % Which session types to use out of above 'session_types'
nSess           = numel(use_sessions);

% what video segments
exp.Video.Segments = {'whiskpadR_wcam','whiskpadL_ecam'}; % 'whiskpadR_ecam','forepawR_ecam', 'whiskpadR_ecam', should be cropped tifs, and suffix indicates camera for timefile.

%% Folder Details (can be imported into workspace) - For Imaging Data

get_local = true;                                                           % READ DATA from local/server? Change to server directory later.
put_local = true;                                                           % SAVE ANALYSIS results in local/server address?

% Folders - local/server addresses - USUALLY CONSTANT for a single user.
% These are 'parent directories' that contain all experiments/animals
exp.Folder.Backup   = 'Y:\Harsha\';                                         % raw data server
exp.Folder.Local    = 'D:\data copy from server\';                       % local copy of data
%'D:\myStuff\data_copy_from_server\'; 

exp.Folder.Analysis = 'Y:\Harsha\Analysed\';                           % analysis directory on server      'X:\Harsha\Golgi in vivo\'
exp.Folder.AnalysisLocal= 'D:\data copy from server\temp analysis\';% local copy of analysis

% << SET >>

% Experiment Details - to change for every experiment. 
% These are 'children' directories WITHIN above parent directories.
exp.Folder.Session  = 'HG 19 Emilie du Chatelet\exp 01\';                     %which experiment/animal/session folder?
exp.Folder.UseRef   = '';                                                   %if patches are to be aligned to REF FROM ELSEWHERE - give << complete path >>
exp.Folder.UseMasks = '';                                                   %if particular masks are to be used - give complete path to .mat

% << IMPORTANT : Analysis will be saved in \\parent dir\children dir\exp.ID >>
exp.ID              = '181128_12_21_49';                                    %Unique ID - Usually represents baseline/ROI/zstack used.
exp.AnimalID        = 1019;                                                 %100N= HG 00N || PIL 1 = Harsha/PIL 2= Fred, 00N = animal number (Harsha's scheme)
exp.Animal          = 'HG 19';
exp.user            = 'Harsha';


% << REMEMBER:   Setup folders for whichever experiments you are using !!!

% Imaging Data (including additional session types)
exp.Folder.baseline = '181128_12_21_49 FunctAcq';                           %baselines can be different for different exp types
exp.Folder.AP       = '181128_12_28_59 FunctAcq';%180409_11_55_20 FunctAcq';   %'180423_12_24_57 FunctAcq'
exp.Folder.AP2      = '181128_12_35_06 FunctAcq';   %'180423_12_29_10 FunctAcq';                    
exp.Folder.baseline2= '';   %'180423_12_38_00 FunctAcq';
exp.Folder.puff_out = '';   %'180423_12_52_57 FunctAcq';
exp.Folder.puff_paw = '';   %'180423_13_21_08 FunctAcq';
exp.Folder.contra_puff= '';
exp.Folder.baseline3= '';   %'180423_13_25_17 FunctAcq';
exp.Folder.extra1   = '';

% Video data
exp.VideoFolder.baseline = '';%180827_11_49_28 VidRec
exp.VideoFolder.AP       = '';%180827_11_58_06 VidRec;%'videos\180423_12_26_36';
exp.VideoFolder.AP2      = '';%'videos\180423_12_17_28';%'videos\180423_12_30_34';
exp.VideoFolder.baseline2= '';%'180827_12_03_53 VidRec';%'videos\180423_12_40_13';
exp.VideoFolder.puff_out = '';%'180827_12_10_31 VidRec';%'videos\180423_12_48_18';
exp.VideoFolder.puff_paw = '';%'videos\180423_13_14_03';%'videos\180423_13_22_39';
exp.VideoFolder.contra_puff ='';
exp.VideoFolder.baseline3= '';%'videos\180423_13_18_39';%'videos\180423_13_27_58';
exp.VideoFolder.extra1   = '';

% ---------------------- Full Paths Made Below Based on Above Settings --------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
if get_local, folder_address = exp.Folder.Local;    else, folder_address = exp.Folder.Backup;         end      % for reading data
if put_local, analysis_address = exp.Folder.AnalysisLocal;  else, analysis_address = exp.Folder.Analysis;   end     % for saving data

for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    exp.(sess)      = [folder_address, exp.Folder.Session, exp.Folder.(sess)];                              % session/exp full paths
end
exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];                          % analysis full path    
mkdir( exp.analysed );

video_split         = true;                                                 %if true, it reads tiff from VideoFolder\split sub-folder, but times from VideoFolder
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    exp.Video.(sess)= [folder_address, exp.Folder.Session, exp.VideoFolder.(sess)];
end

nVid = numel(exp.Video.Segments);
 
exp.session_types = session_types; exp.use_sessions = use_sessions; exp.nSess = nSess;
save( [exp.analysed,'\exp_details.mat'], 'exp', 'session_types', 'nSess', 'use_sessions', '-v7.3');          % save exp details for later use
clear  analysis_address folder_address get_local put_local
%------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------------------------
%% P1 - Basic Pre-processing: (Motion correction, Masking and) Normalisation

exp.params.dff_smooth_scale = 100;
if get_dff    
    run_basic_patch_analysis_multsess();
end
%------- can add option to call patch or points pipeline
%------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------------------------
%% P2 - Analyse behaviour

% Speed
if get_speed
    for sn = 1:nSess
        sess            = session_types{use_sessions(sn)};
        Speed.(sess)    = cell2mat( (get_speed_data([exp.(sess), '\Speed_Data\']))' );
    end

save( [exp.analysed, '\Speed.mat'], 'Speed', '-v7.3' );
disp(['.....................Saved Speed Data'])
end

% Video MI
if get_videoMI
    [ VidMI, Vidtime ] = run_videoMI_multsess(exp);


%Remove duplicates
for sn=1:nSess
    sess = session_types{use_sessions(sn)};
    [~,ia,~] = unique(VidMI.(sess).(Whisk.vidseg)(:,2));
    VidMI.(sess).(Whisk.vidseg) = VidMI.(sess).(Whisk.vidseg)(ia,:);
end
end
% to be added- pupil size, video PCA, behavioural modules, 

%------------------------------------------------------------------------------------------------------------------------------------


%% P3 - Identidying Behavioural Epochs

% << SET >> Parameters for:

% 1. Air Puff times
pufftime=[];
Puff.start  = 2000; %ms                                                     % Puff trigger time every trial
Puff.dur    = 100;  %ms                                                     % Puff duration

% 2. Locomotion Periods using speed
Loco.use_absolute          = true;
Loco.min_period_separation = 500;    % ms
Loco.min_period_size       = 3000;   % ms
Loco.on_threshold          = 0.8;
Loco.off_threshold         = 0.2;
Loco.Avg_window            = 100;    % ms

% 3. Whisking Periods using Motion Index
Whisk.use_absolute          = true;
Whisk.min_period_separation = 500;    % ms
Whisk.min_period_size       = 1000;   % ms
Whisk.on_threshold          = 0.12;
Whisk.off_threshold         = 0.1;
Whisk.Avg_window            = 50;    % ms
Whisk.vidseg                = 'whiskpadR_wcam';

%------------------------------------------------------------------------------------------------------------------------------------
% Compute behaviour analysis
if get_beh_epoch,   get_behaviour_epoch_pipeline;
save( [exp.analysed, '\Beh.mat'], 'pufftime','Loco','Whisk', 'Puff' ,'-v7.3' );
disp('.....................Saved Behaviour Epochs')
end
%------------------------------------------------------------------------------------------------------------------------------------

%% Save distances
ROI.Coordinates = p.POI_coordinates;
ROI.Distance.pw = nan(nROI, nROI);
for roi1=1:nROI
for roi2=roi1+1:nROI
    [ROI.Distance.pw(roi1,roi2), ROI.Distance.pw(roi2,roi1)] = deal( sqrt( (250/512)^2*(ROI.Coordinates.ROIs_X(roi1,1)-ROI.Coordinates.ROIs_X(roi2,1))^2+ (250/512)^2*(ROI.Coordinates.ROIs_Y(roi1,1)-ROI.Coordinates.ROIs_Y(roi2,1))^2 + (ROI.Coordinates.ROIs_Z(roi1,1)-ROI.Coordinates.ROIs_Z(roi2,1))^2));
end
end
save( [exp.analysed,'\ROI.mat'], 'ROI', '-v7.3')

