%% Modified version of Fred and Harsha's motion index code for batch processing
% This version can load the avi and crop ROIs on the fly
% -------------------------------------------------------------------------
% Model : 
%   [motion_indexes, video] = get_MI_from_video(file_path, timestamp, 
%                               rendering, ROI, normalize, dump_data)
%
% -------------------------------------------------------------------------
% Inputs : 
%   file_path_or_data(STR path or analysis_params object or X * Y * T matrix)
%                         - an absolute or relative path pointing to a
%                         multipage image stack or a sequence of eye cam
%                         images. 
%                         - A video as a X * Y * T matrix. dump_data is
%                         automatically set to false in that case
%
%   timestamp(1 * T ARRAY) - Optional - Default 1:Npoints of loaded file
%
%   rendering(BOOL) - Optional - Default true
%                       If true, final result is shown
%
%   ROI(1 * 4 INT or Cell array of {1 * 4} INT) - Optional - Default ''
%                       If specified, the video is cropped and this ROI is
%                       used to calculate the M.
%                       Format is similar to the outut of imrect, which is
%                       [Xstart, Y start, Xsize, Ysize]. You can pass a
%                       cell array of ROIs to get multuple measurment at
%                       once.
%                       If empty, the whole image is used (in case you
%                       cropped the video beforehand)
%
%   normalize(BOOL) - Optional - Default true
%                       If true, motion index is normalized
%
%   dump_data(BOOL) - Optional - Default false
%                       If true, data is stored temporrily on HD to save
%                       memory, but this slow down loading considerably
% -------------------------------------------------------------------------
% Outputs :
%   motion_indexes ({1 * N} CELL ARRAY of [2 * T] ARRAY)
%                       For each ROI input, we get one motion index
%                       measurment. For each cell :
%                       First column is motion index (between 0 and 1).
%                       Second column are timestamps copied from given 
%                       Timestamp vector.
%
%   video (X * Y * T)
%                       The video used for motion index. If you passed the
%                       video as an input, this is the same as
%                       file_path_or_data. If you passed a file path, this
%                       is the corresponding video (UINT8)
% -------------------------------------------------------------------------
% Extra Notes:
% Determine a motion index based on Jelitai et al., Nature communications,
%
% Adapted from Harsha's code. For a more complete version including loading
% from substack and gif, see
% https://github.com/SilverLabUCL/Harsha-old-code/tree/master/Behaviour_analysis
% -------------------------------------------------------------------------
% Author(s):
%   Harsha Gurnani, Frederic Lanore, Antoine Valera
%    
% -------------------------------------------------------------------------
% Revision Date:
%   10-05-2019
%
% See also load_stack
%

function [motion_indexes, video] = get_MI_from_video(file_path_or_data, timestamp, rendering, ROI, normalize, dump_data)
    profile on
    
    %% Open stack
    if isnumeric(file_path_or_data)
        %% If you pass data directly (must be in a X * Y * T format)
        dump_data = false;         
    elseif ~isnumeric(file_path_or_data) && exist(file_path_or_data, 'file')
        %% If you pass a file path and exist, load later
    else
        error('source not identified. Pass a .avi/.tif path, or a [X * Y * T] Matrix \n') 
    end
    if nargin < 3 || isempty(rendering)
        rendering = true;
    end
    if nargin < 5 || isempty(normalize)
        normalize = true;
    end
    if nargin < 6 || isempty(dump_data)
        dump_data = false;
    end
    
    %% Load data if required
    if ~isnumeric(file_path_or_data) && exist(file_path_or_data, 'file')
        %% If you pass a file path 
        fprintf(['please wait... loading videofile :',file_path_or_data,'\n'])
        file_path_or_data = mmread_light(file_path_or_data, dump_data);
        fprintf('video loaded \n')
    end
        
    %% Preparation for ROI collection
    if ~isempty(ROI)
        motion_indexes = cell(1, numel(ROI));
        interbatch_holder = cell(1, numel(ROI));
    else
        motion_indexes = cell(1, 1);
        interbatch_holder = cell(1, 1);
    end
    
    %% Set some useful variables
    if ~dump_data
        %% When video data is available at once (fastest)
        nFrames = size(file_path_or_data, 3);
        file_path_or_data = {file_path_or_data};
        n_src = 1; % When not dumping data, batch size is 1;
        if nargin < 2 || isempty(timestamp)
            timestamp = (1:nFrames)';
        end
    else
        %% When video data is available by batch (memory saving)
        [~, n_src] = size(file_path_or_data,'video_full');
        timestamp = [];
        nFrames = [];
    end
    
    %% Now collect data, batch by batch
    offset = 0;
    
    for batch_idx = 1:n_src  
        %% Setup current batch
        if ~dump_data
            %% Whole video at once
            video = file_path_or_data{batch_idx};
            clear file_path_or_data;
        else
            %% Current batch
            video = file_path_or_data.video_full(1, batch_idx);
            video = video{1};
            nFrames = size(video, 3);
            timestamp = ((1:nFrames) + offset)';
            offset = max(timestamp);
        end
        
        %% Extract MI ROIs to measure, if any
        if nargin >= 4 && ~isempty(ROI)
            if iscell(ROI) % If you passed multiple ROIs as a cell array
                temp = {};
                for el = 1:numel(ROI)
                    roi = round(ROI{el});
                    temp{el} = video(roi(2):roi(2)+roi(4), roi(1):roi(1)+roi(3), :);
                end
                video = temp;
                clear temp
            else % If you passed a single ROI as a 1 * 4 DOUBLE
                video = {video(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3), :)};
            end
        else % If you didn't pass any ROI, then the whole video is your ROI
            video = {video};
        end

        %% Now get MI for each ROI of the current batch
        for MI_idx = 1:numel(video)
            local_data = video{MI_idx};

            %% Preallocate output
            MI = nan(nFrames, 2);  

            %% Compute motion index for first point (for batch_idx > 1)
            if nFrames > 1 && dump_data && batch_idx > 1 
                X  = interbatch_holder{MI_idx};
                X1 = local_data(:,:,1);
                FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));
                MI(1,1) = (sum(sum(FrameDiff.*FrameDiff)));
            end
            %% Compute motion index for the rest
            for i = 2:nFrames
                X  = local_data(:,:,i-1);
                X1 = local_data(:,:,i);
                FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));
                MI(i,1) = (sum(sum(FrameDiff.*FrameDiff)));            
                % FrameDiff = diff(local_data(:,:,i-1:i),1,3).^2;
                % MI(i,1) = sum(FrameDiff(:));
            end
            interbatch_holder{MI_idx} = local_data(:,:,end);

            
            %% Normalize
            if normalize
                m = prctile(MI(:,1),5);
                M = max(MI(:,1));
%                 MI(:,1) = (MI(:,1)-m) ./ (M-m); 
            end

            MI(:,2) = timestamp(1:nFrames,end); %setup time axis
            %MI = MI(2:end,:); 

            %% Stitch to any previous batch
            motion_indexes{MI_idx} = cat(1, motion_indexes{MI_idx}, MI);
        end
    end
    
    %% Render MI's if required
    if rendering 
        figure();hold on
        for el = 1:numel(motion_indexes)
            plot(motion_indexes{el}(:,2),motion_indexes{el}(:,1));hold on;
        end
        %plot(mean(cell2mat(cellfun(@(x) x(:,1), motion_indexes, 'UniformOutput', false)), 2), 'k')
    end
    
    profile off
    profile viewer
end