get_local = true;   % Read Data from local/server? Change to server directory later.
put_local = false;  % Save analysis results in local/server address?


% Folders - local/server addresses - usually constant for a single user
exp.Folder.Backup   = 'W:\Harsha\';
exp.Folder.Analysis = 'X:\Harsha\Golgi in vivo\';
exp.Folder.Local    = 'D:\myStuff\data_copy_from_server\';

% Experiment Details - to change for every experiment
exp.Folder.Session  = 'HG 03 Caroline Herschel\exp 03\';
exp.Folder.Baseline = '180423_12_05_10 FunctAcq';
exp.Folder.AP       = '180423_12_15_54 FunctAcq';%'180319_11_59_31 FunctAcq';
if isempty(exp.Folder.AP), do_AP = false; else, do_AP = true; end

exp.ID              = '180423_12_05_10';
exp.AnimalID        = 1003;    %100N= HG 00N. PIL 1 = Harsha/PIL 2= Fred. 00N = animal number
exp.Animal          = 'HG 03';
exp.user            = 'Harsha';

if get_local, folder_address = exp.Folder.Local;    else, folder_address = exp.DataBackup;         end
if put_local, analysis_address = exp.Folder.Local;  else, analysis_address = exp.Folder.Analysis;   end

exp.baseline        = [folder_address, exp.Folder.Session, exp.Folder.Baseline];
if do_AP
    exp.AP          = [folder_address, exp.Folder.Session, exp.Folder.AP];
end
exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];
suffix.baseline     = 'baseline';
suffix.AP           = 'AP';

mkdir( exp.analysed );
save( [exp.analysed,'\exp_details.mat'], 'exp', '-v7.3');
clear  analysis_address folder_address get_local put_local

%% 1) Basic Pre-processing: Motion correction, masking and Normalisation
run_basic_patch_analysis();
clear suffix 

%% 2) Analyse behaviour

% Speed
Speed.baseline = cell2mat( (get_speed_data([exp.baseline, '\Speed_Data\']))' );
if do_AP
    Speed.AP   = cell2mat( (get_speed_data([exp.AP, '\Speed_Data\']))' );
end

% Plot
nROI = p.maxROI;
cjet = colormap( jet(p.maxROI) );  close all
cj = [cjet(1:2:nROI,:); cjet(2:2:nROI,:)];

plot_exp = 'baseline';
show_patch_and_trace
savefig([exp.analysed,'\',plot_exp,'_norm.fig']);

if do_AP    
    plot_exp = 'AP';
    show_patch_and_trace
    savefig([exp.analysed,'\',plot_exp,'_norm.fig']);
end

disp('.....................Saved Plots')

clear cj plot_exp 

%% 3) Compute correlations

%% 4) determine airpuff responses



