%-------------------------------------------------------------------------%
%                            INSTRUCTIONS                                 %
%-------------------------------------------------------------------------%
% 1. COPY this file in your current directory
% 2. Save it as myAnalysisSettings.m .      DO NOT ERASE OR MODIFY the .template file !!!!!
% 3. Set your parameters in myAnalysisSettings.m. 
% 4. If you want to use the default values, set load_defaults to TRUE. Then, it will NOT read any of your changed parameters, not a single one!
% 5. Once you are happy, run the patch analysis pipeline, make some coffee and let it work for you :)

% ----------------------------------------------------------------------- %



% COMMENT THIS OUT if you want to keep current workspace variables
% Eg if you want to rerun only part of the analysis w
clear variables                                                                                                                                        

init.load_GUI               = true;
init.load_defaults          = false;
init.partial                = false;
init.user_file              = 'setup_harsha.ini';




%% S1 - What Analyses Do You Want

% 1. what measures/recordings?        
exp.get.dff            = true;                                              % (MC/Segmentation for patch data +) Normalisation to get dF/F
exp.get.speed          = true;                                              % Read wheel encoder data
exp.get.videoMI        = false;                                             % Use if video_parts are cropped -> get Motion Index from Videos
exp.get.vidpca         = false;                                             % Dim reduction of motion index/behaviour using PCA
exp.get.ephys          = false;                                             % Simultaneously recorded ephys data (eg spikes)


% 2. what analyses?
exp.get.beh_epoch      = false;                                              % Different Behavioural Epochs
exp.get.corr           = false;                                              % Pairwise and behavioural correlations of dF/F
exp.get.pca            = false;                                              % PCA of dF/F (and background)
exp.get.spike_infer    = false;                                              % Spikes inferred from dF/F (check user_file for algo, SPIKY by default)
exp.get.calibrate      = false;                                             % calibrate spike inference in case of paired ephys and imaging


% 3. what plots?
exp.plot.simple        = true;                                              % Simple plot - (only mean patch,) dff traces and speed
exp.plot.beh           = false;                                              % can add other beh info available - like mi, vid pca etc
exp.plot.puff          = true;                                              % add air puffs
exp.plot.ephys         = false;                                             % Ephys data (eg spikes)
exp.plot.inf_spike     = false;                                             % Inferred spike trains
exp.plot.pw_corr       = false;                                             % Pairwise correlation analysis results
exp.plot.beh_corr      = false;                                              % dff and behaviour correlation analyses results
exp.plot.pca           = false;                                              % PCA results


% 4. what sessions?                                          
exp.setup_type         = 'golgi in vivo';                                   % What kind of experiment - details are in user_file .ini
% Session type names will be used as fields for all analyses and as file suffix:
exp.session_types      =  {'baseline', 'AP', 'AP2', ...                      Eg.        standard exp
                          'baseline2', 'puff_out', 'puff_paw', 'baseline3', ...         control exp
                          'stim_100Hz', 'stim_20Hz', ...                                ephys exp
                          'extra1', 'extra2', 'extra3' };                   % extras are buffers for additional session types/baselines. 
% You can add as many additional session types as you want.  Rename them to relevant descriptor.                                         
exp.use_sessions       = [1 2 ];% 3 4 5 6 7 8 9 10];                       % Which session types to use out of above 'session_types'


% 5. what video segments for MI
exp.Video.Segments     = {'whiskpadR_ecam', 'forepawR_ecam', 'whiskpadR_ecam' };  % should be cropped tifs, and suffix indicates camera for timefile.
exp.Video.TimeFile     = {'ecam',           'ecam',           'ecam'          };  % which timefile, if no _xcam suffix to Video segment. Can be xcam or a text file name
exp.video_split        = true;                                              %if true, it reads tiff from VideoFolder\split sub-folder, but times from VideoFolder
%--------- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
exp.Video.Pupil        = 'pupil_wcam';                                      % format/analysis to follow
exp.Video.PupilTime    = 'wcam';
exp.video_split        = true;                                              %if true, it reads tiff from VideoFolder\split sub-folder, but times from VideoFolder



%% S2 - Folder Details 

% If true, read/write from local address
% else from server - careful about writing to server!!!
exp.get.get_local       = true;                                             % READ DATA from local/server?                  
exp.get.put_local       = true;                                             % SAVE ANALYSIS results in local/server address?

%--------------------------------------------------------------------------
% IMP! Analyses will be saved in    < \parent dir\ children dir\ exp.ID\  >
%--------------------------------------------------------------------------

%------------- META DIRECTORIES -----------------
% can be put in .ini file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ------------- PARENT DIRECTORIES -------------
% These are 'PARENT DIRECTORIES' that contain all experiments/animals
%        - local and server addresses 
%        - USUALLY CONSTANT for a single user.
exp.Folder.read_ini     = false;                                             % Read from user_file .ini file
if ~exp.Folder.read_ini
%for reading
exp.Folder.Backup       = 'W:\Harsha\';                                     % raw data server
exp.Folder.Local        = 'D:\myStuff\data_copy_from_server\';              % local copy of data
% for saving analyses
exp.Folder.Analysis     = 'X:\Harsha\Golgi in vivo\';                       % analysis directory on server
exp.Folder.AnalysisLocal= 'D:\myStuff\data_copy_from_server\Golgi in vivo\';% local copy of analysis
end

% ------------- CHILDREN DIRECTORIES -------------
% These are 'children' directories WITHIN above parent directories.
% Experiment Details - Can be superdirectory for animal/day/FOV
%                    - to change across experiments. 
exp.Folder.Session  = 'HG 13 Emmy Noether\exp 02\';                                      % which experiment/animal/session folder?


% ------------- experiment ID -------------
exp.ID              = '180827_13_05_34';                                    % Unique ID - Usually represents baseline/ROI/zstack used.
                                                                            % Can be any arbitrary string
                                                                            % RECOMMENDED: make it related to session or baseline expt information, or experiment key in your database
% -------------  Optional ------------- 
exp.AnimalID        = 1013;                                                 % UNNN, where U = User Number, NNN = Animal number(Harsha's scheme)
                                                                            % 1003= HG 3 || U=1 -> Harsha; U=2 -> Fred 
exp.Animal          = 'HG 13 Emmy Noether';                            % Animal ID as string -> Name, ID, whatever your naming scheme
exp.user            = 'Harsha';                                             % Who did the experiment
exp.Date            = exp.ID(1:6);                                          % Date of experiment


%----------- SESSION DIRECTORIES -----------------------
% << REMEMBER:   Setup folders for whichever experiments you are using !!!

% IMAGING DATA                                                              (including additional session types)
% Data will be read from < \parent dir\ children dir\ exp.Folder.'session' >
exp.Folder.baseline         = '180827_13_05_34 FunctAcq';                   %baselines can be different for different exp types
exp.Folder.AP               = '180827_12_58_12 FunctAcq';
exp.Folder.AP2              = '';   
exp.Folder.baseline2        = '';  
exp.Folder.puff_out         = '';  
exp.Folder.puff_paw         = '';   
exp.Folder.baseline3        = '';   
exp.Folder.stim_100Hz       = '';
exp.Folder.stim_20Hz        = '';
exp.Folder.extra1           = '';
% ADD exp.Folder.NewSessionName for any additional session type 'NewSessionName' you use.


% VIDEO DATA
% Data will be read from < \parent dir\ children dir\ exp.VideoFolder.'session' >
exp.VideoFolder.baseline    = '';
exp.VideoFolder.AP          = '';
exp.VideoFolder.AP2         = '';
exp.VideoFolder.baseline2   = '';
exp.VideoFolder.puff_out    = '';
exp.VideoFolder.puff_paw    = '';
exp.VideoFolder.baseline3   = '';
exp.VideoFolder.extra1      = '';


% EPHYS DATA
% Data will be read from < \parent dir\ children dir\ exp.EphysFolder.'session' >
exp.EphysFolder.baseline   = '';
exp.EphysFolder.stim_100Hz = '';
exp.EphysFolder.stim_20Hz  = '';





%% S3 - Other Files for Processing
exp.Use.Ref           = '';                                                 % if patches are to be aligned to REF FROM ELSEWHERE - give << complete path >> to a .mat file 
                                                                            % Needs to be a workspace with all relevant variables, with at least a Structure called 'MC' with field 'Ref'
exp.Use.Masks         = '';                                                 % if particular masks are to be used - give complete path to .mat file                                                                            
exp.Use.Grouping      = '';                                                 % If multiple ROI belong to same cell - can be cell array of integers (for only somatic patches) OR
                                                                            % Struct array with fields soma, dend, axon, background, upper, lower, other

                                                                            



%% S4 - Pre Processing 

% raw data extraction
Basic.ROI_type            = 'patch'; % or 'points', 'plane'                 % What kind of data did you record?

Basic.ROIs                = 0;                                              % which ROI to keep for analysis;    0 = All, else give array with ROI numbers
Basic.Groups              = 0;                                              % Which groups to keep for analysis; 0 = all
Basic.trials              = 0;                                              % Which trials to keep for analysis; 0 = all

Basic.continuous_time     = true;                                           % Single continuous time - computed based on cycle lengths* number of timepoints
                                                                            % inter-trial time can be added if speed file is available.
                                                                            % If false, each trial is considered a repetition
% motion correction (for patches and plane timeseries)
Basic.MC                  = true;                                           % If true, corr-based motion correction will be applied to planes/patches before any segmentation/normalisation
Basic.MC_select_ROI       = true;                                           % If true, user input is required at the time of execution to select ROI used for MC. For patches, that is patch#, for planes, user can draw rectangles.
                                                                            % Can be set to an array (for patch#) or 5-column matrix (each row is a rectangle - [plane #, x1, xn, y1, yn]) instead of boolean value
                                                                            % If false, all patches are used for MC - can be memory intensive!!
Basic.MC_max_disp         = 8;                                              % Max displacement allowed in any MC-frame; If optimal disp> max_disp, F = NaN in that frame                                                                            
Basic.MC_show             = false;                                          % Show movie (random frames) of Movement-corrected data

% segmentation/averaging (for patches and plane timeseries)
Basic.patch_avg           = false;                                          % If true, mask(s) is applied to patch, else each pixel is normalised independently.
Basic.patch_multiple      = false;                                          % If true, multiple masks are kept for a single patch(plane)

Basic.Mask_MinSize        = 10;                                             % Min number of pixels in accepted mask
Basic.get_background      = true;                                           % Identify background_per_patch.
Basic.ROI_background      = [];                                             % Array of ROI numbers that were put as backgrounds

% normalisation and smoothing
Basic.norm_baseline       = 'percentile';    % or 'interval'                % How to compute F0?
Basic.base_interval       = [];%ms                                          % Time interval considered as baseline, and its time-averaged F = F0
Basic.base_percentile     = 10;                                             % percentile to fet F0 per trace
Basic.base_invert         = false;                                          % If true, events are assumed to be negative, and hence normalisation is changed accordingly. 
                                                                            % Top base_percentile rather than lowest base_percentile is considered as baseline
Basic.baseline_pertrial   = false;                                          % if true, each trial has its own baseline

Basic.smooth_scale        = 100;%ms                                         % Filter size for smoothing normalised traces (boxcar_filter)
Basic.smooth_causal       = false;                                          % If true, filter is causal i.e [t-n : t] rather than [t-n/2 : t+n/2]

% for ephys data
Basic.ephys_spikes        = 'cell-attached';                                % If not empty, spikes are extracted from ephys data depending on the recording method specified





%% S5 - Identidying Behavioural Epochs/ Other Events

% Events can be externally triggered or determined from data.
% This is useful for event-triggered or event-averaged analyses
% Check Event Type for more details on fields.

Event.types                       = {'Puff', 'Loco','Whisk'};  %,'Spike'    % Add Event Names to be used as fieldnames and file suffixes, and then add details below for analyses


% 1. Air Puff times
Event.Puff.type                   = 'ext';                                  % 'ext' - external with times defined, or 'beh' - find behaviour epoch, or 'spike' - use ephys or inference data
% params events with defined times
Event.Puff.times                  = [];                                     % If not empty, 2-column matrix with start and stop times of all events (third column can be repetition if in trial time)
                                                                            % Can also be a .mat file address
Event.Puff.start                  = 2000; %ms                               % If 'times' are empty, event trigger time every trial ( scalar or array with time for each trial)
Event.Puff.dur                    = 100;  %ms                               % Puff duration                 ( scalar or array with time for each trial)
Event.Puff.func_eval              = @airpuff_analysis;

% 2. Locomotion Periods using speed
Event.Loco.type                  = 'beh';
% params for events to identify
Event.Loco.use_absolute          = true;                                    % use absolute value of encoder data to find 'locomotion'?
Event.Loco.min_period_separation = 500;    % ms                             % min separation between 2 identified beh periods, or else they are merged (including time between)
Event.Loco.min_period_size       = 3000;   % ms                             % min length of identified beh period
Event.Loco.on_threshold          = 0.8;                                     % threshold of moving average for start of beh period
Event.Loco.off_threshold         = 0.2;                                     % threshold of moving average for end of beh period
Event.Loco.Avg_window            = 100;    % ms                             % window for moving average
Event.Loco.func_eval             = @run_speed_multsess;

% 3. Whisking Periods using Motion Index
Event.Whisk.type                  = 'beh';
Event.Whisk.use_absolute          = true;
Event.Whisk.min_period_separation = 500;    % ms
Event.Whisk.min_period_size       = 1000;   % ms
Event.Whisk.on_threshold          = 0.12;
Event.Whisk.off_threshold         = 0.1;
Event.Whisk.Avg_window            = 50;    % ms
Event.Whisk.vidseg                = 'whiskpadR_ecam';                       % which whisker video to use for whisking periods (in case of multiple videos or tracking)
Event.Whisk.func_eval             = @run_videoMI_multsess;


%% S6 - Activity Correlations

% 1. PAIRWISE correlations for dF/F:  < Corr.pw.params >
% ----------------------------------
% Computed for total trace per session, and also during different
% behavioural periods specified
Corr.pw.ds_rate             = 50;           % Hz                            % Rate at which all data is sampled to compute pairwise correlations. Time bin of correlation is effectively 1/ds_rate
Corr.pw.maxlag              = 2000;         % ms                            % Maximum lag for cross correlation
Corr.pw.nshuffle            = 500;                                          % No. of shuffles to compute null distribution of correlation coefficients

Corr.pw.event               = {'Loco', 'Whisk', 'Puff'};                    % Event Epochs to compute PW in
Corr.pw.Loco.extratime      = 500;          % ms                            % Extra time around each identified beh period (locomotion here) to use for computing correlation for each beh period
                                                                            % Can be scalar (same extra time before and after beh period) OR array: [before after] extra time
Corr.pw.Whisk.extratime     = 500;          % ms                            % Extra time around whisking periods
Corr.pw.Puff.extratime      = [200 1000];   % ms                            % Extratime around air puff period


% 2. Correlation WITH BEHAVIOUR :                 < Corr.beh.params >
% -------------------------------
% correlation between dF/F and behavioural trace
% Computed for total trace and identified event epochs
Corr.beh.behaviour          = {'Speed', 'WhiskMI'}; %   , 'FR'}             % What traces to correlated dF/F with?                   
Corr.beh.ds_rate            = 10;           % Hz
Corr.beh.maxlag             = 6000;         % ms
Corr.beh.nshuffle           = 500;

Corr.beh.event              = {'Loco', 'Whisk_noLoco', 'Whisk_noPuff', 'Puff'};
Corr.beh.Loco.extratime     = 500;          %ms
Corr.beh.Whisk_noLoco.extratime = 200;      %ms
Corr.beh.Whisk_noPuff.extratime = 200;      %ms
Corr.beh.Puff.extratime     = [200 1000];   %ms





%% S7 - Activity PCA

PCA.numComp                  = 10;                                          % Max number of population dF/F PCs to keep. 
PCA.minVarExp                = [];                                          % If numComp is empty, keep as many PCs as need to explain 'minVarExp' percentage of variance           

% Behavioural correlations ( of top numComp PCs with behavioural traces
PC_Corr.params.ds_rate       = 20;
PC_Corr.params.maxlag        = 2000;
PC_Corr.params.nshuffle      = 500;
% Events are picked up from Corr structure



%% S8 - Plotting

% 1. Traces
% ----------------
Plot.BasicPlot.two_col          = true;
Plot.BasicPlot.scale            = 30;                                       % 0 if it should be determined automatically
Plot.BasicPlot.showPatch        = true;
Plot.BasicPlot.graylim          = [0 500];                                  % 0 if it should be determined automatically - for gray patch mean image
Plot.BasicPlot.sortDepth        = 'ascend';                                 % Ascend means deepest neurons are lowest traces. 'decend' is reverse. If empty. then it doesn't sort.

Plot.BasicPlot.LineWidth        = 1.5;
Plot.BasicPlot.ColorMap         = 'jet';                                    % Can be a matrix
Plot.BasicPlot.EventColor.Loco  = 'k';
Plot.BasicPlot.EventColor.Whisk = 'g';
Plot.BasicPlot.EventColor.Puff  = 'm';


% 2. PCA results
% ----------------
Plot.PCA_Plot.expVar            = {'baseline', 'all' };                     % which sessions? for ExpVar v/s NumComp
Plot.PCA_Plot.trace             = {'baseline', 'all' };                     % Top n principal components of which session
Plot.PCA_Plot.nComp             = 5;                                        % How many components on the PCA trace plot? If 0, use as many as saved

Plot.PCA_Plot.ColorMap          = 'jet';
Plot.PCA_Plot.LineWidth         = 1.5;


% 3. PW Correlation
% ------------------
%   Plot PW somatic corr of type PW_Plot.type
%   total in session 1, 2, ... , and for event1.. in session 1, etc .
Plot.PW_Plot.type               = {'zerolag', 'max'};                       % What kind of PW correlation matrix to plot?
Plot.PW_Plot.sess               = { { 'baseline', 'Whisk', 'Loco'}, ...     { {'session1', 'event1','event2'}, {'session2', 'event1', ..}, {'session3', 'event1',..}, ... }
                                    { 'AP', 'Whisk', 'Loco', 'Puff'}...     % If event list is empty, then it only plots total PW corr. 
                                  };                                        % Can also use 'all' as session/event name

% 3. Behavioural Correlation
% --------------------------                          
Plot.BCor_Plot.corrstyle        = 'trace';  % or 'imgsc'                    % corr vs lag plot
Plot.BCor_Plot.hist_type        = { 'zerolag', 'max' };
Plot.BCor_Plot.dfftype          = { 'soma', 'pc' };

Plot.BCor_Plot.sess             = { { 'baseline', 'whisk', 'loco', 'pupil' }, ... { {'session1', 'beh1','beh2'}, {'session2', 'beh1', ..}, {'session3', 'beh1',..}, ... }
                                    { 'AP', 'whisk', 'loco', 'pupil'}...   Plot somatic corr with beh1&2 in session 1, with beh1.. in session 2, with beh1... in session3 ...
                                  };                                                   

