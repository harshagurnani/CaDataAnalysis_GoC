myAnalysisSettings;

%% Set up parameters/folders
% Change myAnalysisSettings.m or load exp into workspace, or continue and
% use GUI
if ~exist('exp','var'),     exp = get_params_gui('exp');    end
if ~isfield(exp, exp.session_types{exp.use_sessions(1)})
    exp = make_complete_address( exp, exp.get.get_local, exp.get.put_local);
end
save_all_params( exp, Basic, Event, Corr, PCA , Plot );

%------------------------------------------------------------------------------------------------------------------------------------
%% P1-A     Basic Pre-processing: (Motion correction, Masking and) Normalisation

% essential fields :
%   In exp:     'get.dff', 'session_types', 'use_sessions', session folder addresses, analysis address
%   In Basic:   'ROI_type' (rest will be set to default)
if exp.get.dff
    [Basic, Seg, Norm] = run_basic_dff_multsess( exp, Basic );
end

%------------------------------------------------------------------------------------------------------------------------------------
%% P1-B 	Extract behaviour traces

% essential fields :
% In exp:     'get.speed', session types and folder addresses,
if exp.get.speed
   Speed = run_speed_multsess( exp );
end

% Video MI
% In exp:     'get.videoMI', session_types, video folder addresses, 'video_split'
if exp.get.videoMI
    [ VidMI, Vidtime ] = run_videoMI_multsess(exp);
end

% to be added- pupil size, 
if ~isempty( exp.Video.Pupil)
    [ Pupil ] = run_pupilAnalysis_multsess(exp); 
end

%% P1-C Extract Task Events (to be added for tasks later)


%% Save distances

test_session = exp.session_types{exp.use_sessions(1)};
ROI.Coordinates = Norm.(test_session).p.POI_coordinates; 
nROI = size(ROI.Coordinates.ROIs_X,1);
ROI.Distance.pw = nan(nROI, nROI);
for roi1=1:nROI
for roi2=roi1+1:nROI
    [ROI.Distance.pw(roi1,roi2), ROI.Distance.pw(roi2,roi1)] = deal( sqrt( (250/512)^2*(ROI.Coordinates.ROIs_X(roi1,1)-ROI.Coordinates.ROIs_X(roi2,1))^2+ (250/512)^2*(ROI.Coordinates.ROIs_Y(roi1,1)-ROI.Coordinates.ROIs_Y(roi2,1))^2 + (ROI.Coordinates.ROIs_Z(roi1,1)-ROI.Coordinates.ROIs_Z(roi2,1))^2));
end
end
save( [exp.analysed,'\ROI.mat'], 'ROI', '-v7.3')
clear test_session nROI 
disp('.....................Saved ROI Coordinates')

%------------------------------------------------------------------------------------------------------------------------------------
%%                       HELPER FUNCTIONS                               
%------------------------------------------------------------------------------------------------------------------------------------
function exp = make_complete_address(exp, get_local, put_local)
% Full Paths made on settings in exp structure

    if get_local, folder_address = exp.Folder.Local;    else, folder_address = exp.DataBackup;         end      % for reading data
    if put_local, analysis_address = exp.Folder.AnalysisLocal;  else, analysis_address = exp.Folder.Analysis;   end     % for saving data

    exp.nSess = length(exp.use_sessions);
    for sn = 1:exp.nSess
        sess            = exp.session_types{exp.use_sessions(sn)};
        exp.(sess)      = [folder_address, exp.Folder.Session, exp.Folder.(sess)];                              % session/exp full paths
        exp.Video.(sess)= [folder_address, exp.Folder.Session, exp.VideoFolder.(sess)];
    end
    exp.analysed        = [analysis_address, exp.Folder.Session, 'analysed\', exp.ID];                          % analysis full path    
    mkdir( exp.analysed );

    exp.nVid = numel(exp.Video.Segments);
    
    save( [exp.analysed,'\exp_details.mat'], 'exp', '-v7.3'); % save exp details for later use

end


function save_all_params( exp, Basic, Event, Corr, PCA , Plot )
% save original settings

    save([exp.analysed, '\Params_exp.mat'], 'exp', '-v7.3');
    save([exp.analysed, '\Params_basic.mat'], 'Basic', '-v7.3');
    save([exp.analysed, '\Params_event.mat'], 'Event', '-v7.3');
    save([exp.analysed, '\Params_corr.mat'], 'Corr', '-v7.3');
    save([exp.analysed, '\Params_pca.mat'], 'PCA', '-v7.3');
    save([exp.analysed, '\Params_plot.mat'], 'Plot', '-v7.3');

end