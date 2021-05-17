clear variables
partial = false;

%% S1 - What Analyses Do You Want

% what measures?        
exp.get.dff            = true;                                              % (MC/Segmentation for patch data +) Normalisation to get dF/F
exp.get.speed          = true;                                              % Read wheel encoder data
exp.get.videoMI        = true;                                              % Use if video_parts are cropped -> get Motion Index from Videos

exp.get.vidpca         = false;                                             % Dim reduction of motion index/behaviour using PCA

% what analyses
exp.get.beh_epoch      = true;                                              % Different Behavioural Epochs
exp.get.corr           = true;                                              % Pairwise and behavioural correlations of dF/F
exp.get.pca            = true;                                              % PCA of dF/F (and background)

% what plots?
exp.get.plot_simple    = true;                                              % Simple plot - (only mean patch,) dff traces and speed
exp.get.plot_beh       = true;                                              % can add other beh info available - like mi, vid pca etc
exp.get.puff           = true;                                              % add air puffs

% what sessions?
exp.session_types      = {'baseline', 'AP', 'AP2', ...                      %these will be used as suffix identifiers
                          'baseline2', 'puff_out', 'puff_paw', 'baseline3', ...
                          'extra1', 'extra2', 'extra3' };                   % extras are buffers for additional session types - or different baselines. You can add as many additional session types as you want.
                                                                            % Rename them to relevant descriptor.                                                                           
exp.use_sessions       = [1  2 ];% 3 4 5 6 7];                              % Which session types to use out of above 'session_types'

% what video segments
exp.Video.Segments     = {'whiskpadR_ecam', 'forepawR_ecam', 'whiskpadR_ecam' };% should be cropped tifs, and suffix indicates camera for timefile.
exp.video_split        = true;                                              %if true, it reads tiff from VideoFolder\split sub-folder, but times from VideoFolder


%% S2 - Folder Details (can be imported into workspace) - For Imaging Data

%------------- META DIRECTORIES -----------------

% If true, read/write from local address
% else from server - careful about writing to server!!!
exp.get.get_local       = true;                                             % READ DATA from local/server?                  
exp.get.put_local       = true;                                             % SAVE ANALYSIS results in local/server address?

%           < parent dir >
% Folders - local and server addresses 
%         - USUALLY CONSTANT for a single user.
% These are 'PARENT DIRECTORIES' that contain all experiments/animals
exp.Folder.Backup       = 'W:\Harsha\';                                     % raw data server
exp.Folder.Local        = 'D:\myStuff\data_copy_from_server\';              % local copy of data

exp.Folder.Analysis     = 'X:\Harsha\Golgi in vivo\';                       % analysis directory on server
exp.Folder.AnalysisLocal= 'D:\myStuff\data_copy_from_server\Golgi in vivo\';% local copy of analysis

%           < children dir >
% Experiment Details - Can be superdirectory for animal/day/FOV
%                    - to change across experiments. 
% These are 'children' directories WITHIN above parent directories.
exp.Folder.Session  = 'HG 03\exp 02\';                                      % which experiment/animal/session folder?

%           < experiment ID >
% IMPORTANT: Analysis will be saved in < \parent dir\ children dir\ exp.ID >
exp.ID              = '180409_11_45_58';                                    % Unique ID - Usually represents baseline/ROI/zstack used.
                                                                            % Can be any arbitrary string
                                                                            % RECOMMENDED: make it related to session or baseline expt information, or experiment key in your database
%       < Optional: >
exp.AnimalID        = 1003;                                                 % UNNN, where U = User Number, NNN = Animal number(Harsha's scheme)
                                                                            % 1003= HG 3 || U=1 -> Harsha; U=2 -> Fred 
exp.Animal          = 'HG 03 Caroline Herschel';                            % Animal ID as string -> Name, ID, whatever your naming scheme
exp.user            = 'Harsha';                                             % Who did the experiment
exp.Date            = exp.ID(1:6);                                          % Date of experiment


%----------- SESSION DIRECTORIES -----------------------
% << REMEMBER:   Setup folders for whichever experiments you are using !!!

% IMAGING DATA                                                              (including additional session types)
% Data will be read from < \parent dir\ children dir\ exp.Folder.'session' >
exp.Folder.baseline         = '180409_11_45_58 FunctAcq';                   %baselines can be different for different exp types
exp.Folder.AP               = '180409_11_55_20 FunctAcq';
exp.Folder.AP2              = '';   
exp.Folder.baseline2        = '';  
exp.Folder.puff_out         = '';  
exp.Folder.puff_paw         = '';   
exp.Folder.baseline3        = '';   
exp.Folder.extra1           = '';
% ADD exp.Folder.NewSessionName for any additional session type 'NewSessionName' you use.


% VIDEO DATA
% Data will be read from < \parent dir\ children dir\ exp.VideoFolder.'session' >
exp.VideoFolder.baseline    = '180409_11_52_30 VidRec\180409_11_52_30 VidRec';
exp.VideoFolder.AP          = '180409_11_56_43 VidRec\180409_11_56_43 VidRec';
exp.VideoFolder.AP2         = '';
exp.VideoFolder.baseline2   = '';
exp.VideoFolder.puff_out    = '';
exp.VideoFolder.puff_paw    = '';
exp.VideoFolder.baseline3   = '';
exp.VideoFolder.extra1      = '';
exp.VideoFolder.extra2      = '';

%% S3 - Pre Processing Settings
exp.Folder.UseRef           = '';                                           % if patches are to be aligned to REF FROM ELSEWHERE - give << complete path >> to a .mat file 
                                                                            % Needs to be a workspace with all relevant variables, with at least a Structure called 'MC' with field 'Ref'
exp.Folder.UseMasks         = '';                                           % if particular masks are to be used - give complete path to .mat file
                                                                            


%% S4 - Identidying Behavioural Epochs

% Parameters for:

% 1. Air Puff times
pufftime                    = [];
Puff.start                  = 2000; %ms                                     % Puff trigger time every trial ( scalar or array with time for each trial)
Puff.dur                    = 100;  %ms                                     % Puff duration                 ( scalar or array with time for each trial)

% 2. Locomotion Periods using speed
Loco.use_absolute          = true;                                          % use absolute value of encoder data to find 'locomotion'?
Loco.min_period_separation = 500;    % ms                                   % min separation between 2 identified beh periods, or else they are merged (including time between)
Loco.min_period_size       = 3000;   % ms                                   % min length of identified beh period
Loco.on_threshold          = 0.8;                                           % threshold of moving average for start of beh period
Loco.off_threshold         = 0.2;                                           % threshold of moving average for end of beh period
Loco.Avg_window            = 100;    % ms                                   % window for moving average

% 3. Whisking Periods using Motion Index
Whisk.use_absolute          = true;
Whisk.min_period_separation = 500;    % ms
Whisk.min_period_size       = 1000;   % ms
Whisk.on_threshold          = 0.12;
Whisk.off_threshold         = 0.1;
Whisk.Avg_window            = 50;    % ms
Whisk.vidseg                = 'whiskpadR_ecam';                             % which whisker video to use for whisking periods (in case of multiple videos or tracking)


%% S5 - Activity Correlations

% PAIRWISE correlations for dF/F:  < Corr.pw.params >
% Computed for total trace per session, and also during different
% behavioural periods
Corr.pw.params.ds_rate             = 50;           % Hz                     % Rate at which all data is sampled to compute pairwise correlations. 
                                                                            % Time bin of correlation is effectively 1/ds_rate
Corr.pw.params.maxlag              = 2000;         % ms                     % Maximum lag to compute cross correlation for
Corr.pw.params.nshuffle            = 500;                                   % No. of shuffles to compute null distribution of correlation coefficients
Corr.pw.params.loco.extratime      = 500;          % ms                     % Extra time around each identified beh period (locomotion here) to use for computing correlation for each beh period
                                                                            % Can be scalar (same extra time before and after beh period) OR array: [before after] extra time
Corr.pw.params.whisk.extratime     = 500;          % ms                     % Extra time around whisking periods
Corr.pw.params.puff.extratime      = [200 1000];   % ms                     % Extratime around air puff period

% WITH BEHAVIOUR :                 < Corr.beh.params >
% correlation between dF/F and behavioural trace
% Computed for total trace ONLY (for now)
Corr.beh.params.ds_rate            = 10;           % Hz
Corr.beh.params.maxlag             = 6000;         % ms
Corr.beh.params.nshuffle           = 500;

%% S6 - Activity PCA

PCA.numComp                  = 10;                                           % Max number of population dF/F PCs to keep 

% Behavioural correlations ( of top numComp PCs with behavioural traces
PC_Corr.params.ds_rate       = 20;
PC_Corr.params.maxlag        = 2000;
PC_Corr.params.nshuffle      = 500;


%% S7 - Plotting

% Traces
PlotParams.two_col          = true;
PlotParams.scale            = 30;
PlotParams.graylim          = [0 500];

% PCA results
PlotParams.PCA_expVar       = {'baseline', 'all' };                             % which sessions?
PlotParams.PCA_trace        = {'baseline', 'all' };                             % Top n components of which session
PlotParams.nComp            = 5;                                                % How many components?

% Correlation plots
PlotParams.corrstyle        = 'trace';  % or 'imgsc'
% of individual ROI
PlotParams.soma_sess        = { { 'baseline', 'whisk', 'loco', 'pupil'      }, ... { {'session1', 'beh1','beh2#}, {'session2', 'beh1', ..}, {'session3', 'beh1',..}, ... }
                                { 'AP', 'whisk', 'ap_whisk', 'loco', 'pupil'}...   Plot somatic corr with beh1,2 in session 1, with beh1.. in session 2, with beh1... in session3 ...
                            };                                                   
% of top PC components
PlotParams.PCA_sess         = { { 'baseline', 'whisk', 'loco', 'pupil'      }, ... { {'session1', 'beh1','beh2#}, {'session2', 'beh1', ..}, {'session3', 'beh1',..}, ... }
                                { 'AP', 'whisk', 'ap_whisk', 'loco', 'pupil'}...   Plot somatic corr with beh1,2 in session 1, with beh1.. in session 2, with beh1... in session3 ...
                            };                                                   

