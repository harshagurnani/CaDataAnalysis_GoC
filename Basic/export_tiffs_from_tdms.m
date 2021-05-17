%% Set up path and folders to analyse

% LOCATION OF MICROSCOPE CONTROLLER
controller_path  = 'D:\Matlab repos\microscope_controller';          

% EXPERIMENT FOLDERS
data_dir =  'E:\Golgi in vivo imaging\Raw Data\Win_dual_MF_L'; %Folder with experiment folders to export

%        Uncomment 1 if you want ALL experiments in data_dir to be exported
% OR     Uncomment 2 if you want a SINGLE experiment in data_dir to be exported

%  1 - MULTIPLE EXPERIMENTS -----------------
allexp = dir([data_dir,'\*FunctAcq*']);     % Look for all LABVIEW folders

% %   % 2 -  SINGLE EXPERIMENT -------------------
% exp_folder = '200103_15_25_41 FunctAcq';    % Actual folder
% allexp = dir([data_dir,'\*',exp_folder]);


% OTHER PARAMETERS
nTrialsLoad = 30;       % Number of trials to process
trial_batch_size = 3;   % reduce it if 'OUT OF MEMORY ERROR' --> will be fixed later 
                        % trials are loaded in batches of this size only for exporting to .mat
                        %.tiff will still have all repeats 

do_save_tiff = false;    % Switch to false if you only want to export .mat files

do_transfer_mat = true; % Also transfer the folder to 'Matlab-like-export-folder'

% number of files should be read from header files
% transfer mat to directory noby itself

%% ------------------------  PROCESSING     ----------------------------%%
addpath(genpath( controller_path ));     
cd( controller_path );
numexp = numel(allexp);
% 
for fn = 1:numexp
    fprintf(' ..................................................................\n .................    Exporting Experiment # %d of %d\n', fn, numexp) 
    labview_exp = [data_dir,'\',allexp(fn).name];
    C = load_ini_file([labview_exp, '\Experiment Header.ini']);
    nTrialsAll = read_ini_value(C,'number of trials',read_ini_value(C,'Number of trials',[], '[FUNCTIONAL IMAGING]'));
    nTrials = min(nTrialsLoad, nTrialsAll)

    % Save trial_batch_size # of trials at a time - until the saving is fixed
    for jj=1:ceil(nTrials/trial_batch_size)
        repeats  = (1:trial_batch_size) + (jj-1)*trial_batch_size; repeats = repeats(repeats<=nTrials);
        p=analysis_params('source',labview_exp, 'data_type', 'raw', 'repeats', repeats, ...         % 5 trials at a time
                          'signal_channel', 2, 'linearize', true, ...                               % only green channel for tiffs
                          'branch_rendering',false, 'rendering', false, 'branch_analysis',false);   % turn viewer off
        load_TDMS( p, 'source',labview_exp, 'RAW', true, 'rendering', false);                       %create .mat exports
    end
    
    matlab_name_format = strrep( allexp(fn).name(8:15), '_', '-');
    
    matlab_exp = [data_dir, '\', matlab_name_format];    
    p = analysis_params('source',matlab_exp, 'repeats', 0, 'signal_channel', 2, 'linearize', true, ...% only green channel for tiffs
                          'branch_rendering',false, 'rendering', false, 'branch_analysis',false, 'smoothing', [0 0]);   % turn viewer off
    % Write tiffs
%     if do_save_tiff
%         for jj=1:params.n_rois
%             repeats  = (1:trial_batch_size) + (jj-1)*trial_batch_size; repeats = repeats(repeats<=nTrials);
%             p = analysis_params('source',matlab_exp, 'repeats', 0, 'ROIs', jj, 'signal_channel', 2, 'linearize', true, ...% only green channel for tiffs
%                           'branch_rendering',false, 'rendering', false, 'branch_analysis',false, 'smoothing', [0 0]);   % turn viewer off
%             save_tiff( p );
%         end
%     end
        
    if do_transfer_mat
        newdir = [labview_exp, '\Matlab-like-export-folder\'];
        movefile( [matlab_exp], [newdir ,'\', strrep( allexp(fn).name(8:15), '_', '-')])
    end
    
end

%% ------------------------ HELPER FUNCTION ----------------------------%%
function save_tiff( p )

    fprintf('...exporting data as .tif file ... please wait\n');
    export_folder = [p.data_folder, '\Matlab_Export\'];
    mkdir(export_folder);
    
    [vol,~,data] = load_experiment(p);
    data = permute(data,[1,2,3,4,6,7,5]); %put the channel in the end
    % data = [x y z time roi trial channel ]
    h2 = waitbar(0,'Saving Tiffs ...');
    vol.signal_channel=1;
    for roi =  vol.ROIs
        roi_str = ['_ROI_' sprintf('%02d',roi)];
        
        d = data(:,:,:,:,roi,:,vol.signal_channel);
        d = reshape(d,size(d,1),size(d,2),[],1);
        waitbar(roi \ numel(vol.ROIs),h2,['Saving Tif for ROI # ', sprintf('%d ', roi)])
        save_stack(d,[export_folder, 'Matlab_Export',roi_str,'.tif']);
    end
    fprintf('...tif file saved\n');
    close(h2)
    
end