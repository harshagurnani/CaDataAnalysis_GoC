clear variables
get_local = true;   % Read Data from local/server? Change to server directory later.
put_local = false;  % Save analysis results in local/server address?

% extras are buffers for additional session types - or different baselines. Rename them to relevant descriptor.
session_types   = {'baseline', 'AP', 'puff_out', 'puff_paw', 'AP2', 'baseline2', 'baseline3', 'extra1,', 'extra2' };    %these will be used as suffix identifiers 
use_sessions    = [1 ];%2 3 4 5 6 7];    %Which session types out of session_types
nSess           = numel(use_sessions);

exp.Folder.Session  = 'FL 90\exp 01\';                  %which experiment/animal/session folder?
exp.ID              = '180316_16_00_05';    %Unique ID - Usually represents baseline/ROI marking used.
exp.Folder.Analysis = 'X:\Harsha\Golgi in vivo\';
exp.Folder.Local    = 'D:\myStuff\data_copy_from_server\';  

% set up directories for reading/writing data
folder_address      = exp.Folder.Local;
analysis_address    = exp.Folder.Analysis;
exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];

load( [exp.analysed,'\exp_details.mat'], 'exp')
load( [exp.analysed, '\masked_norm_test.mat'], 'Norm', 'p')
load( [exp.analysed, '\masked_raw.mat'], 'Seg')

%% Video directories
exp.VideoFolder.baseline = 'videos\180316_16_05_53';       %baselines can be different for different exp types
exp.VideoFolder.AP       = 'videos\180423_12_26_36';%'180319_11_59_31 FunctAcq';
exp.VideoFolder.AP2      = 'videos\180423_12_30_34';
exp.VideoFolder.baseline2= 'videos\180423_12_40_13';
exp.VideoFolder.puff_out = 'videos\180423_12_48_18';%'180319_11_59_31 FunctAcq';
exp.VideoFolder.puff_paw = 'videos\180423_13_22_39';%'180319_11_59_31 FunctAcq';
exp.VideoFolder.baseline3= 'videos\180423_13_27_58';
exp.VideoFolder.extra1   = '';
exp.VideoFolder.extra2   = '';

video_split              = true;    %if true, it reads tiff from videofolder\split subfolder, but times from videofolder
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    exp.Video.(sess)      = [folder_address, exp.Folder.Session, exp.VideoFolder.(sess)];
end

exp.Video.Segments = {'forepawR_ecam', 'whiskpadR_ecam', 'whiskpadR_wcam'};
nVid = numel(exp.Video.Segments);

%% Get Motion Index
run_videoMI_multsess;

%% Plot and Save
% Plot

nROI = numel(Seg.masks);
cjet = colormap( jet(nROI) );  close all
cj = [cjet(1:2:nROI,:); cjet(2:2:nROI,:)];
% Speed
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    Speed.(sess)    = cell2mat( (get_speed_data([exp.(sess), '\Speed_Data\']))' );
end

% add_puff = false;
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};

    
    % get puff times:
    if ~isempty(strfind(sess, 'AP')) || ~isempty(strfind(sess,'puff') )
        add_puff = true;
        tmp = Norm.(sess).t;    nTr = 10;   nTm = size(tmp,1)/nTr; 
        tmp = reshape(tmp,[nTm,nTr,nROI]);
        puff.start  = 2000;%ms
        puff.dur    = 100;  %ms
        pufftime.(sess) = puff.start + arrayfun( @(trial) min(tmp(1,trial,:)), 1:nTr)-min(tmp(:));  %ms
        
    else
        add_puff = false;
    end
    show_patch_and_videotrace
    savefig([exp.analysed,'\',sess,'_vidMI.fig']);

end

% disp('.....................Saved Plots')

clear cj sess 
