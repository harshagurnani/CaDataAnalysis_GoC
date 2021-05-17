function [fileList] = saveImagesFromVideo(videoFile, saveFolder, varargin)

%Saves single .png grayscale images from a videoFile into a pre-existing  
%saveFolder. The starting and ending time as well as the timeInterval 
%between the images can be specified. Additionally, if only a section of 
%the image should be saved a bounding box can be supplied and the image can 
%be scaled to decrease size. Returns and stores a fileList into the 
%saveFolder with the path of every image and the CurrentTime used to 
%extract the image.
%   
%INPUT 
%   videoFile = File from which the videos should be extracted. 
%   saveFolder = Folder into which the data will be saved. 
%
%OPTIONAL INPUT - pass as a structure with the appropriate fields
%   startTime = Time in s (any precision) at which image
%       extraction starts. Default = 0.
%   endTime = Time to which the images are extracted. If this is larger
%       than the lenth of the video it will be capped to the length. 
%       Default = length of the video. 
%   timeInterval = Interval (in s, any precision) between 
%       frames, if this is smaller than the video speed it will default
%       to the video speed. Default = video speed.
%   boundingBox = y,x - section of the video that should be saved. Size of
%       boudingBox has the format [yStart xStart height width].
%   imageScaling = value between 0 and 1 denoting the scaling factor for
%       imresize.
%
%OUPUT
%   Single .png grayscale images will be saved into the specified folder,
%   the naming is frame1.png, frame2.png etc. 
%   fileList = list of the full file names (including path) where the 
%       images were stored and the CurrentTime used to extract the image. 
%       fileList is also saved into the saveFolder as 
%       'frameNameAndTime.mat'

%% Validate Input
assert(exist(videoFile, 'file') == 2, 'Video file not found')
if ~(saveFolder(length(saveFolder)) == '\')
    saveFolder = strcat(saveFolder, '\');
end
assert(exist(saveFolder, 'dir') == 7, 'Target folder not found')
assert(~(exist(strcat(saveFolder, 'frame1.png'), 'file') == 2), ...
    'File with the name frame1.png already exists in the saveFolder')

%% Extract video data and add optional parameters

camera.Object = VideoReader(videoFile);
startTime = (1/camera.Object.FrameRate)/2;
endTime = camera.Object.Duration;
timeInterval = 1/camera.Object.FrameRate;
imageScaling = 1;

if ~isempty(varargin)
    optionalInput = varargin{1}; 
    if isfield(optionalInput, 'startTime')
        validateattributes(optionalInput.startTime, {'numeric'}, ...
            {'nonnegative'}, 'saveImagesFromVideo', 'startTime')
        startTime = optionalInput.startTime;
        assert(startTime < endTime, 'startTime does not fit video length')
    end
    if isfield(optionalInput, 'endTime')
        assert(optionalInput.endTime <= endTime, ...
            'endTime larger than the length of the video')
        endTime = optionalInput.endTime;
        assert(endTime >= startTime, 'endTime does not fit with startTime')
    end
    if isfield(optionalInput, 'timeInterval')
        validateattributes(optionalInput.timeInterval, {'numeric'}, ...
            {'positive'}, 'saveImagesFromVideo', 'timeInterval')
        assert(optionalInput.timeInterval >= timeInterval, ...
            'Chosen timeInterval is smaller than recording interval')
        timeInterval = optionalInput.timeInterval;
    end
    if startTime == 0 
        startTime = timeInterval/2;
        % Else the first image might be doubled (frame1 & frame2)
    end
    if isfield(optionalInput, 'boundingBox')
        validateattributes(optionalInput.boundingBox, {'numeric'}, ...
            {'positive','size',[1 4]}, 'saveImagesFromVideo','boundingBox')
        boundingBox = optionalInput.boundingBox;
    end
    if isfield(optionalInput, 'imageScaling')
        validateattributes(optionalInput.imageScaling, {'numeric'}, ...
            {'positive', '<=', 1}, 'saveImagesFromVideo', 'imageScaling')
        imageScaling = optionalInput.imageScaling;
    end
end

%% Read video frames
frameCount = 0;
for iTime = startTime:timeInterval:endTime
    camera.Object.CurrentTime = iTime;
    videoFrame = readFrame(camera.Object);
    if size(videoFrame,3)>1
        videoFrame = rgb2gray(videoFrame);
    end  
    
    %% Adjust size if boundingBox given
    if exist('boundingBox', 'var')
        newYStart = boundingBox(1,1); 
        newXStart = boundingBox(1,2);
        newHeight = boundingBox(1,3); 
        newWidth = boundingBox(1,4);
        videoFrame = videoFrame(newYStart : newYStart+newHeight-1, ...
            newXStart : newXStart+newWidth-1);
    end
    if imageScaling < 1
        videoFrame = imresize(videoFrame, imageScaling);
    end
    
    %% Write single image
    frameCount = frameCount + 1;
    imwrite(videoFrame, strcat(saveFolder, 'frame', ...
        int2str(frameCount), '.png'))
    fileList(frameCount).name = strcat('frame',int2str(frameCount),'.png');
    fileList(frameCount).Time = iTime;
end
save(strcat(saveFolder, 'frameNameAndTime.mat'), 'fileList')
