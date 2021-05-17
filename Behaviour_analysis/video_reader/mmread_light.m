function video_full = mmread_light(filename, dump_data)
    if nargin < 2 || isempty(dump_data)
        dump_data = false;
    end

    trySeeking = true;
    matlabCommand = '';
    
    %% Quick query to get some useful values & Array preallocation
    vidObj = VideoReader(filename);
    numFrames = ceil(vidObj.FrameRate*vidObj.Duration);
    batchsize = 500;
    batches = [1:batchsize:numFrames,numFrames+1]; %% FFGrab('doCapture') uses a lot of memory  if we load everythin at once, so we do smaller batches;
    idx_range = 1:3:(vidObj.Width * vidObj.Height * 3); 
    video_full = cell(1, numel(batches(1:end-1)));    
    if dump_data
        delete('temp_video.mat');
        save('temp_video.mat' ,'video_full', '-v7.3','-nocompression')
        file = matfile('temp_video.mat','Writable',true);
    end
    
    %% Data collection, for each batch
    for idx = 1:numel(batches(1:end-1))
        fprintf([num2str(100*idx/numel(batches(1:end-1))),'%% done\n']);
        frames = batches(idx):(batches(idx+1)-1);
        FFGrab('build',filename,'',double(false),double(true),double(trySeeking));
        FFGrab('setFrames',frames);
        FFGrab('setMatlabCommand',matlabCommand);
        FFGrab('doCapture');

        [nrVideoStreams, nrAudioStreams] = FFGrab('getCaptureInfo');

        
        %% Loop through getting all of the video data from each stream

        for i=1:nrVideoStreams
            [width, height, rate, nrFramesCaptured, nrFramesTotal, totalDuration] = FFGrab('getVideoInfo',i-1);
            if (nrFramesTotal > 0 && any(frames > nrFramesTotal))
                warning('mmread:general',['Frame(s) ' num2str(frames(frames>nrFramesTotal)) ' exceed the number of frames in the movie.']);
            end

            %% Memory preallocation
            video = zeros(height * width * 3, nrFramesCaptured, 'uint8');
            
            %% Get each frame
            for f = 1:nrFramesCaptured
                [video(:,f), ~] = FFGrab('getVideoFrame',i-1,f-1);
                % the data ordering is wrong for matlab images, so permute it
            end
            
            %% keep one channel only and reshape
            if ~isempty(video)
                video = permute(reshape(video(idx_range(idx_range <= size(video, 1)),:), width, height, nrFramesCaptured),[2 1 3]); % WARNING -- Disabled channel 2 and 3 here
                if dump_data 
                    file.video_full(1, idx) = {video};
                end
            end
        end
        FFGrab('cleanUp'); 
        if ~dump_data
            video_full{idx} = video;
        end
    end
    
    if ~dump_data
        %% Contatenate along time axis
        video_full = cat(3, video_full{:});
    else
        video_full = matfile('temp_video.mat','Writable',false);
    end
end
