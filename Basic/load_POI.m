function [params, raw_data, timedata, norm_data] = load_POI(varargin)
%% Load ROIs automatically using files in path <fullpath> - with Gui for selecting load parameters
% ----------------------
% Usage : 
%
% [params,raw_data, timedata, norm_data] = load_POI()                         --- Pop-up appears for inputing options
%           OR
% [params,raw_data, timedata, norm_data] = load_POI('Name1', 'Value1', .... ) --- No pop-up (except for folder unless that option is given). Name-value
%                                                                           arguments are parsed to modify those arguments and rest are set to default.
% ----------------------
%
% Available Options (as Name-Value pairs):
% NAME                      TYPE                DEFAULT VALUE           DETAILS
%--------------------------------------------------------------------------------
% 'source'                  (STR)               ''                      Full experiment folder path, otherwise pop-up appears
% 'ROIs'                    0 or INT array      0                       ROI ids to load and analyse ( 0 = All ROI)
% 'ROI_type'                (STR)               'poi'                   'poi' or 'patch'
% 'add_itt'                 (BOOL)              True                    Add inter-trial time (ITT) (from speed file)
% 'avg_patch'               (BOOL)              True                    Average raw F of whole patch
% 'trials'                  0 or INT array      0                       Trial numbers to load and analyse (0 = All trials)
% 'dwell_time'              (FLOAT)             ''                      Dwell time in ms. Defaults to 0.0001 for ROI_type = 'patch'
% 'channel'                 (STR)               'green'                 'red' or 'green'
% 'tdms2int'                (BOOL)              True                    Multiply rawF by dwell_time/(5 ns) -> LabView data averages rather than
%                                                                       adds data collected in 5 ns bins.
% 'normalize'               (BOOL)              True                    Normalise as DFF?
% 'pop_up'                  (BOOL)              True if no option give, otherwise false. Pop-up for parameters.
%
%
% If 'normalize' == true, following options are additionally useful.
% NAME                      TYPE                DEFAULT VALUE           DETAILS
%--------------------------------------------------------------------------------
% 'norm_method'             (STR)               'percentile'            Normalization method - What is the baselinne? Options: Median, Percentile
% 'base_percentile'         INT of FLOAT        10                      Percentile value of baseline for norm_method= 'Percentile' 
% 'baseline_per_trial'      (BOOL)              false                   Compute new baseline for each trial?
% 'smooth_type'             (STR)               'Moving average'        Type of smoothing. (Only one optiona currently available)                                          
% 'smooth_scale'            (INT of FLOAT)      70                      Time/extent of square window for smoothing (in ms)

% ------------------- Not functional options yet
% {'save_norm', (BOOL) },                   If true, save normalized data
% {'speed_data',(STR)},                     Path of speed log data
% {'flatten',(BOOL)},                       If true, Volmes are flattened
% {'colormap',(cmap object)},               matlab colormap
% %%%%{'save',(BOOL)},                      if false, concat/averaged stacks are not saved
% 
% 
% OUTPUT:
%-------------------------
% params     - Structure with import parameters     
%
% [raw_data] - a (Timepoints x Trial# x ROIs#) 
%              Each 2D array (:,:,ii) is one ROI, with columns as different trials.
% [timedata] - Same size as rawdata. Correct timestamps for rawdata
% [normdata] - DFF: [Timepoints x Trial# x ROIs#) normalized result 

[ getpath, params, popup_on, get_coord ] = parse_options(nargin, varargin, 0);

if getpath
params.exp_path = uigetdir(pwd, 'Choose experiment directory...');
params.exp_path = strrep(params.exp_path, '\', '/');
end


%% Import parameters
if popup_on
    prompt = {'ROI_type (poi or patch)','ROIs (0 = all)','Trials (0 = all)','Dwell time (ms)',...
           'Fill time (ms)','green or red', 'Average values across patch?', 'Multiply patch tdms data by dwell time bins?',  'Add time between trials?', 'Normalize and smoothen?'};
    dlg_title = 'Importing parameters';
    num_lines = 1;

    %Default params
    defaultans = {params.ROI_type, num2str(params.ROIs), num2str(params.trials), num2str(params.dwell_time), num2str(params.fill_time), params.signal_channel, num2str(params.average_patch), num2str(params.tdms_to_int), num2str(params.add_intertime), num2str(params.normalize) };
    ans_ROIs = inputdlg(prompt,dlg_title,num_lines,defaultans);

    params.ROI_type = ans_ROIs{1,1};
    params.ROIs = sort(str2num(ans_ROIs{2,1}));
    params.trials = sort(str2num(ans_ROIs{3,1}));
    params.dwell_time = str2num(ans_ROIs{4,1}); % in ms.
    params.fill_time = str2num(ans_ROIs{5,1}); % in ms.
    params.signal_channel = ans_ROIs{6,1};
    params.average_patch =  strcmp( ans_ROIs{7,1}, 'true' ) || strcmp( ans_ROIs{7,1}, '1' );
    params.tdms_to_int = strcmp( ans_ROIs{8,1}, 'true' ) || strcmp( ans_ROIs{8,1}, '1' );
    params.add_intertime =strcmp( ans_ROIs{9,1},'true') || strcmp( ans_ROIs{9,1},'1');
    params.normalize =strcmp( ans_ROIs{10,1},'true') || strcmp( ans_ROIs{10,1},'1'); 
end

h = load_ini_file([params.exp_path,'/Experiment Header.ini']);
params = read_params_from_header( h, params );


%Coordinates - LabView format - ROI.dat file
POI_coordinates = importdata([params.exp_path, '/ROI.dat']);
switch params.ROI_type
    case 'poi'
        POI_coordinates.ROIs_X = POI_coordinates.data(:,5);
        POI_coordinates.ROIs_Y = POI_coordinates.data(:,6);
        POI_coordinates.ROIs_Z = POI_coordinates.data(:,7);
        
        if exist([params.exp_path, '/Functional_Data/'],'dir')
            params.data_path = [params.exp_path, '/Functional_Data/'];
            params.acq_software = 'Labview'; params.import_folder = 'Labview';
        elseif exist([params.exp_path, '/Matlab-like-export-folder'], 'dir')
            subf = ls([params.exp_path, '/Matlab-like-export-folder/']);
            subf = subf(end,:);
            params.data_path = [params.exp_path, '/Matlab-like-export-folder/', subf, '/'];    
            params.acq_software = 'Labview'; params.import_folder = 'Matlab-new';
        else
            % for recordings with matlab.
        end
    case 'patch'
        POI_coordinates.ROIs_X = [ POI_coordinates.data(:,5), POI_coordinates.data(:,8) ];
        POI_coordinates.ROIs_Y = [ POI_coordinates.data(:,6), POI_coordinates.data(:,9) ];
        POI_coordinates.ROIs_Z = [ POI_coordinates.data(:,7), POI_coordinates.data(:,10) ];
        
        %From Antoine's exporter. To be made robust to further changes to
        %name of export directory structure.
        if exist([params.exp_path, '/Exported_Experiment_1/'], 'dir')
            params.data_path = [params.exp_path, '/Exported_Experiment_1/ThisShouldBeTime/'];  
            params.acq_software = 'Labview'; params.import_folder = 'Matlab-old';
        elseif exist([params.exp_path, '/Matlab-like-export-folder'], 'dir')
            subf = ls([params.exp_path, '/Matlab-like-export-folder/']);
            subf = subf(end,1:8);
            if isfolder([params.exp_path, '/Matlab-like-export-folder/',subf])
                params.data_path = [params.exp_path, '/Matlab-like-export-folder/', subf, '/'];  
            else
                params.data_path = [params.exp_path, '/Matlab-like-export-folder/'];
            end
            params.acq_software = 'Labview'; params.import_folder = 'Matlab-new';
        else
            %Read TDMS
        end
        
end

if get_coord
    params.POI_coordinates = POI_coordinates;
end

%Total ROI scanned
try
cycletimes = importdata([params.exp_path '/Single cycle relative times.txt']); % values are in microseconds.
catch

cycletimes = [];
end
params.maxROI = size(POI_coordinates.ROIs_X, 1);


if params.ROIs == 0
    %case: all ROIs
    params.nROIs = params.maxROI;
    params.ROIs = 1:params.nROIs;
else
    params.nROIs = length(params.ROIs);
end


switch params.ROI_type
    case 'poi'
        if strcmp(params.import_folder,'Labview')
            numfiles = length(dir( [params.data_path, '*POI*.dat'] ))/2;    %2 files (red/green channels) per POI x trial
        else
            numfiles = length(dir( [params.data_path, '*Scan*.mat'] ));  %Both channels in 1 file
        end
    case 'patch'
         numfiles = length(dir( [params.data_path, '*Scan*.mat'] ));  %Both channels in 1 file
end
params.maxTrials = numfiles/params.maxROI;
if params.trials == 0
    %case: all trials
    params.nTrials = params.maxTrials;
    params.trials = 1:params.nTrials;
else
    params.nTrials = length(params.trials);
end

%% Smoothing parameters
if params.normalize && popup_on
        prompt = {'Normalization method','Baseline percentile','baseline per trial?','Smoothing type', 'Temporal smoothing scale (ms)'};
        dlg_title = 'Normalization and Smoothing parameters';
        num_lines = 1;

        %Default params
        defaultans = {params.norm_method, num2str(params.base_percentile), num2str(params.baseline_per_trial), params.smooth_type, num2str(params.smooth_scale) };
        ans_smoothing = inputdlg(prompt,dlg_title,num_lines,defaultans);

        params.norm_method = ans_smoothing{1,1};
        params.base_percentile = str2num(ans_smoothing{2,1});
        params.baseline_per_trial = strcmp( ans_smoothing{3,1}, 'true') || strcmp( ans_smoothing{3,1}, '1');
        params.smooth_type = ans_smoothing{4,1}; 
        params.smooth_scale = str2num(ans_smoothing{5,1}); % in ms
end

%% Get timestamps and raw data
disp('Importing raw data...')
% Convert patch data from mean to sum
if  strcmp(params.ROI_type, 'patch') && params.tdms_to_int
    if isempty( params.dwell_time) || params.dwell_time==0
        disp( ' No value for dwell time - using 100 ns by default ' );
        params.dwell_time = 0.0001;
    end
    params.mult_factor = floor(params.dwell_time/ 5e-6 );
    disp( sprintf('Multiplying raw data by %d to change averages to sums ...', params.mult_factor) )     %#ok<DSPS>
end
tic
[raw_data, params] = get_rawdata(params);
toc



disp('Assigning timestamps...')
tic
 [timedata, params] = get_timestamps2(params, cycletimes);
toc

%% Normalize to get df/F0, and smoothen

if params.normalize
    disp('Normalising to dF/F...')
    tic
    if strcmp(params.ROI_type, 'patch') && ~params.average_patch        
        ts = size(raw_data);
        temp =  reshape( permute( raw_data, [3,4,5,1,2] ) , [ts(3), ts(4), ts(5)*ts(1)*ts(2)] );
        time_temp = reshape( permute( timedata, [3,4,5,1,2] ) , [ts(3), ts(4), ts(5)*ts(1)*ts(2)] );
        [norm_data,F0] = normalize_and_smooth(temp, time_temp, params );
        clear temp time_temp
        norm_data = permute( reshape( norm_data, [ts(3), ts(4), ts(5), ts(1), ts(2)]), [4, 5, 1, 2, 3] );
        params.F0 = permute( reshape( F0, [ts(4), ts(5), ts(1), ts(2)]), [3,4,1,2] );
    else
        [norm_data, params.F0] = normalize_and_smooth(raw_data, timedata, params);
    end
    
    toc
else
    norm_data = [];
end
end


function [ getpath, params, popup_on, get_coord ] = parse_options( nargs, varargs, copts )

    getpath=true;
    params.ROI_type = 'poi';
    params.ROIs = 0;
    params.trials = 0;
    params.dwell_time = 0.0001; %num2str(0.0001);
    params.fill_time = 0.0245;
    params.signal_channel = 'green';
    params.average_patch = true;
    params.tdms_to_int = true;
    params.add_intertime=true;
    params.normalize=true;

    params.norm_method = 'percentile';
    params.base_percentile = 10;
    params.baseline_per_trial = false;
    params.smooth_type = 'Moving average';
    params.smooth_scale = 100;

    popup_on=true;
    get_coord = true;%false;
    
if nargs > copts  
    popup_on = false;
    start_read_var = copts+1;
    while start_read_var <= nargs
       vn = start_read_var-copts;
       if ischar( varargs{vn} )
           if strcmp(varargs{vn}, 'source' )
               params.exp_path = varargs{vn+1};
               getpath=false;
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'ROIs' )
               params.ROIs = varargs{vn+1};
               start_read_var = start_read_var + 2;     
           elseif strcmp(varargs{vn}, 'ROI_type' )
               params.ROI_type = varargs{vn+1};
               start_read_var = start_read_var + 2;     
           elseif strcmp(varargs{vn}, 'add_itt' )
               params.add_intertime = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'avg_patch' )
               params.average_patch = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'trials' )
               params.trials = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'dwell_time' )
               params.dwell_time = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'channel' )
               params.signal_channel = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'tdms2int' )
               params.tdms_to_int = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn},  'normalize' )
               params.normalize = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'pop_up' )
               popup_on = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'norm_method' )
               params.norm_method = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'base_percentile' )
               params.base_percentile = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'baseline_per_trial' )
               params.baseline_per_trial = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'smooth_type' )
               params.smooth_type = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'smooth_scale' )
               params.smooth_scale = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'get_coord')
               get_coord = varargs{vn+1};
               start_read_var = start_read_var + 2;
           else 
                error('Argument name %s is not an optional argument', varargs{vn})
           end
       end
    end

end
end



function params = read_params_from_header( h, params )

    new_software = ~isnan(read_ini_value( h, 'Software Version', NaN, '[LOGIN]'));
    
    if new_software
        
        % Acq details:
        
        % Imaging Mode
        
        
        % Number of trials
        if params.trials == 0
            nTrials = read_ini_value( h, 'Number of trials', 0, '[FUNCTIONAL IMAGING]');
            params.trials = 1:nTrials;
        end
        % Number of timepoints
        params.nTimepoints = read_ini_value( h, 'Number of cycles', 0, '[FUNCTIONAL IMAGING]');
    
        
        
        % FOV
        params.FOV = read_ini_value( h, 'field of view', 250, '[FUNCTIONAL IMAGING]');  %um
        params.nFramePixels = read_ini_value( h, 'Frame Size', 512, '[FUNCTIONAL IMAGING]');   
        params.xy_res = params.FOV/params.nFramePixels;
        sx = read_ini_value( h, 'Stage x position', 0, '[FUNCTIONAL IMAGING]');  %um
        sy = read_ini_value( h, 'Stage y position', 0, '[FUNCTIONAL IMAGING]');  %um
        sz = read_ini_value( h, 'Stage z position', 0, '[FUNCTIONAL IMAGING]');  %um
        params.stagePos = [sx sy sz];
        
        
        % Roi details
        if strcmpi(params.ROI_type, 'patch')
            % Number of lines/pixels
            params.nPatchLines = read_ini_value( h, 'Rectangle ROI height', 1, '[FUNCTIONAL IMAGING]');
            params.nLinePixels = read_ini_value( h, 'Rectangle ROI length', 1, '[FUNCTIONAL IMAGING]');
            % Number of ROIs
            if params.ROIs == 0
                nTotalLines = read_ini_value( h, 'Number of miniscans', 1, '[FUNCTIONAL IMAGING]');
                params.ROIs = 1: (nTotalLines/params.nPatchLines);
            end
            % Dwell time
            params.dwell_time = read_ini_value( h, 'pixel dwell time (us)', 1, '[FUNCTIONAL IMAGING]')/1000;    %us->ms

        elseif strcmpi(params.ROI_type, 'poi')
            % Dwell time
            params.dwell_time = read_ini_value( h, 'poi dwell time', 1, '[FUNCTIONAL IMAGING]')/1000;   %us->ms
            % Number of ROIs
            if params.ROIs == 0
                params.ROIs = 1: read_ini_value( h, 'Number of pois', 1, '[FUNCTIONAL IMAGING]');
            end
        end


    else
        
    end
    
end
