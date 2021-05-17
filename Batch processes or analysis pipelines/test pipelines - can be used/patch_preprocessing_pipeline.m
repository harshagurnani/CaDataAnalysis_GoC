clear variables
get_local = true;   % Read Data from local/server? Change to server directory later.
put_local = false;  % Save analysis results in local/server address?

% extras are buffers for additional session types - or different baselines. Rename them to relevant descriptor.
session_types   = {'baseline', 'AP', 'puff_out', 'puff_paw', 'AP2', 'baseline2', 'baseline3' };    %these will be used as suffix identifiers 
use_sessions    = [1 2 ];% 3 4 5 6 7];    %Which session types out of session_types
nSess           = numel(use_sessions);

% Folders - local/server addresses - usually constant for a single user
exp.Folder.Backup   = 'W:\Harsha\';
exp.Folder.Analysis = 'X:\Harsha\Golgi in vivo\';
exp.Folder.Local    = 'D:\myStuff\data_copy_from_server\';  

% Experiment Details - to change for every experiment
exp.Folder.Session  = 'FL 90\exp 02\';                  %which experiment/animal/session folder?
exp.Folder.UseRef   = '';                               %if patches are to be aligned to ref from elsewhere - give complete address
exp.Folder.UseMasks = '';                               %if particular masks are to be used - give complete mat address

exp.Folder.baseline = '180427_11_43_48 FunctAcq';       %baselines can be different for different exp types
exp.Folder.AP       = '';%'180423_12_24_57 FunctAcq';%'180319_11_59_31 FunctAcq';
exp.Folder.puff_out = '';%'180423_12_52_57 FunctAcq';%'180319_11_59_31 FunctAcq';
exp.Folder.puff_paw = '';%'180423_13_21_08 FunctAcq';%'180319_11_59_31 FunctAcq';
exp.Folder.AP2      = '';%'180423_12_29_10 FunctAcq';
exp.Folder.baseline2= '';%'180423_12_38_00 FunctAcq';
exp.Folder.baseline3= '';%'180423_13_25_17 FunctAcq';

exp.ID              = '180427_11_43_48';    %Unique ID - Usually represents baseline/ROI marking used.
exp.AnimalID        = 2090;                 %100N= HG 00N || PIL 1 = Harsha/PIL 2= Fred, 00N = animal number
exp.Animal          = 'FL 90';
exp.user            = 'Harsha';

% set up directories for reading/writing data
if get_local, folder_address = exp.Folder.Local;    else, folder_address = exp.DataBackup;         end
if put_local, analysis_address = exp.Folder.Local;  else, analysis_address = exp.Folder.Analysis;   end

for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    exp.(sess)      = [folder_address, exp.Folder.Session, exp.Folder.(sess)];
end
exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];
    
mkdir( exp.analysed );
save( [exp.analysed,'\exp_details.mat'], 'exp', '-v7.3');
clear  analysis_address folder_address get_local put_local

%% 1) Basic Pre-processing: Motion correction, masking and Normalisation
run_basic_patch_analysis_multsess();


%% 2) Analyse behaviour
% Speed
for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    Speed.(sess)    = cell2mat( (get_speed_data([exp.(sess), '\Speed_Data\']))' );
end

% Plot
nROI = p.maxROI;
cjet = colormap( jet(p.maxROI) );  close all
cj = [cjet(1:2:nROI,:); cjet(2:2:nROI,:)];

for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    plot_exp        = sess;

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
    
    show_patch_and_trace
    savefig([exp.analysed,'\',plot_exp,'_norm.fig']);

end

disp('.....................Saved Plots')

clear cj plot_exp 

%% 3) Compute correlations

%% 4) determine airpuff responses



