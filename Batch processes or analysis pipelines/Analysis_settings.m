%% load parameters
% SET ALL in myAnalysisSettings file OR Load Defaults
% For partial analysis: Write over, but get/set settings should not be written over. Only
% associated with file.
try
    myAnalysisSettings
catch
    load_default = true; partial = false;
    disp('No file called myAnalysisSettings.m foumd. Loading defaults...')
end
if load_default,    load_analysis_defaults;     end


% Params Structure - for all necessary/ ubiquitously used parameters;
Params.nSess           = numel(use_sessions);
Params.Sessions        = arrayfun( @(n) exp.session_types{ n}, exp.use_sessions, 'UniformOutput', false );

%% Set complete Paths 
if get_local, folder_address = exp.Folder.Local;    else, folder_address = exp.DataBackup;         end      % for reading data
if put_local, analysis_address = exp.Folder.AnalysisLocal;  else, analysis_address = exp.Folder.Analysis;   end     % for saving data

% for imaging
for sn = 1:Params.nSess
    sess            = Params.Sessions(sn);
    exp.(sess)      = [folder_address, exp.Folder.Session, exp.Folder.(sess)];                              % session/exp full paths
end
exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];                          % analysis full path    
mkdir( exp.analysed );

% for videos
for sn = 1:Params.nSess
    sess            = Params.Sessions(sn);
    exp.Video.(sess)= [folder_address, exp.Folder.Session, exp.VideoFolder.(sess)];
end
Params.nVid = numel(exp.Video.Segments);


%% save analysis settings
if partial
    exptemp = exp;
    load([exp.analysed,'\exp_details.mat'], 'exp');
    tempget = exp.get;  exp = exptemp;
    exp.get = tempget;
    clear tempget exptemp
end

save( [exp.analysed,'\exp_details.mat'], 'exp', 'Params', '-v7.3');         % save exp details for later use
clear  analysis_address folder_address get_local put_local

try
   copyfile( 'myAnalysisSettings.m', [exp.analysed,'\']);
catch
end