function videoOfTracedCenters(tracedCenters, tracedImages, saveFolder, ...
    saveFile, varargin)
%Saves a uncompressed AVI video file of the supplied tracedImages with the
%corresponding tracedCenters indicated by a 3x3 pixel red dot. Can add a
%second set of traced centers with optinal input.
%
%INPUT
%   tracedCenters = a n-by-2 matrix of the y,x coordinates of the traced
%       feature in n tracedImages.
%   tracedImages = the uint8 grayscale images corresponding to the traced
%       Centers. Can besupplied as a folder name in which the images
%       "frame1.png", "frame2.png" ... will be treated as belonging to
%       tracedCenters(1,:), tracedCenters(2,:)...
%       Alternatively can be supplied as a struct with the field
%       'imageFilename'. The images will be assumed to be in the same order
%       as the tracedCenters.
%   saveFolder = folder into which the video file is saved, full path name
%   saveFile = name under which the file is saved, needs to end in '.avi'
%
%OPTINAL INPUT - pass as a structure with the appropriate fields
%   tracedCenters2 = a n-by-2 matrix of the y,x coordinates of a second set
%       of features traced in the same images as tracedCenters
%       tracedImages. Default = no second set.
%   frameRate = integer of the frame rate that the video is saved with.
%       Default = 10.
%   subset = specifies the start and end frame [start end] of a subset if
%       not all tracedCenters and corresponding tracedImages should be part
%       of the video. Default = no subset. Even if a subset is supplied it
%       is expected that tracedCenters(1,:) corresponds to 'frame1.png' or
%       the first image in the list.
%   useFirstImage = specifies whether the first image (frame1.png) should
%       be used or not. If tracing was done frame1.png might not have a
%       corrsesponding traceCenter point. 1 = use the first image, 0 = omit
%       the first image. Default = 1.
%   displayFrameNumber = specifies whether the frame number should be
%       displayed in the upper left corner of every frame. 1 = display, 0 =
%       don't display. Default = 0. Line to change font size or location
%       is indicated in the code.
%
%OUTPUT
%   Video file stored into the saveFolder with the saveFile name.

%% Validate Inputs
validateattributes(tracedCenters, {'numeric'}, ...
    {'size', [NaN, 2], 'positive'}, 'videoOfTracedCenters','tracedCenters')
assert(~isempty(strfind(saveFile, '.avi')) && strfind(saveFile, '.avi')...
    == length(saveFile)-3, 'saveFile does not end with .avi')
assert(isdir(saveFolder), 'Target folder not found')
if ~(saveFolder(length(saveFolder)) == '\')
    saveFolder = strcat(saveFolder, '\');
end
assert(~exist(strcat(saveFolder, saveFile), 'file'), ...
    'Target file already exists')
if ischar(tracedImages) && isdir(tracedImages)
    if ~(tracedImages(length(tracedImages)) == '\')
        tracedImages = strcat(tracedImages, '\');
    end
    frameList = dir(strcat(tracedImages, 'frame*.png'));
    for iFrame = 1:length(frameList)
        frameNames(iFrame).imageFilename = strcat(tracedImages, 'frame',...
            int2str(iFrame), '.png');
    end
elseif isstruct(tracedImages)
    frameNames = struct('imageFilename', {tracedImages.imageFilename});
else
    error('tracedImages not found')
end
assert(~isempty(frameNames), 'tracedImages not found')


%% Get frameRate,tracedCenters2,displayFrameNumber & pos. remove first img
if ~isempty(varargin)
    optionalInput = varargin{1};
    frameRate = 10;
    if isfield(optionalInput, 'frameRate')
        frameRate = optionalInput.frameRate;
    end
    if isfield(optionalInput, 'tracedCenters2')
        tracedCenters2 = optionalInput.tracedCenters2;
        assert(size(tracedCenters2) == size(tracedCenters), ...
            'the two center traces do not have the same size')
    end
    displayFrameNumber = 0;
    if isfield(optionalInput, 'displayFrameNumber')
        displayFrameNumber = optionalInput.displayFrameNumber;
    end
    if isfield(optionalInput, 'useFirstImage') && ...
            optionalInput.useFirstImage == 0
        frameNames(1) = [];
    end
    assert(length(frameNames) == length(tracedCenters), ...
        'There are not the same number of images and traced centers')
    
    %% Implement subset values if supplied
    if isfield(optionalInput, 'subset')
        subset = optionalInput.subset;
        assert(subset(1,1)<subset(1,2) && subset(1,1)>0 && ...
            subset(1,2)<length(tracedCenters), ...
            'subset selection does not match the number of traced centers')
        tracedCenters = tracedCenters(subset(1,1):subset(1,2), :);
        frameNames(subset(1,2)+1:length(frameNames)) = [];
        if subset(1,1) > 1
            frameNames(1:subset(1,1)-1) = [];
        end
        if exist('tracedCenters2', 'var')
            tracedCenters2 = tracedCenters2(subset(1,1):subset(1,2), :);
        end
    end
end

%% Read in image files and tag traced positions
imageFrame = imread(frameNames(1).imageFilename);
tracedCenters = round(tracedCenters);
videoFrames = zeros(size(imageFrame,1), size(imageFrame,2) , 3, ...
    length(frameNames), 'uint8');
for iFrame = 1:length(frameNames)
    imageFrame = imread(frameNames(iFrame).imageFilename);
    videoFrames(:,:,1,iFrame) = imageFrame;
    videoFrames(:,:,2,iFrame) = imageFrame;
    videoFrames(:,:,3,iFrame) = imageFrame;
    if sum(tracedCenters(iFrame,:))>0
        frameY = tracedCenters(iFrame,1);
        frameX = tracedCenters(iFrame,2);
        videoFrames(frameY-1:frameY+1, frameX-1:frameX+1, 1, iFrame) = 255;
        videoFrames(frameY-1:frameY+1, frameX-1:frameX+1, 2, iFrame) = 0;
        videoFrames(frameY-1:frameY+1, frameX-1:frameX+1, 3, iFrame) = 0;
    end
end
if exist('tracedCenters2', 'var')
    tracedCenters2 = round(tracedCenters2);
    for iFrame = 1:length(frameNames)
        if sum(tracedCenters2(iFrame,:))>0
            frameY = tracedCenters2(iFrame,1);
            frameX = tracedCenters2(iFrame,2);
            videoFrames(frameY-1:frameY+1, frameX-1:frameX+1, 2, ...
                iFrame) = 255;
            videoFrames(frameY-1:frameY+1, frameX-1:frameX+1, 3, ...
                iFrame) = 0;
        end
    end
end

%% Insert frame number if displayFrameNumber = 1
if displayFrameNumber == 1
    for iFrame = 1:length(frameNames)
        % Change location and size of displayed numbers here
        videoFrames(:,:,:,iFrame) = insertText...
            (videoFrames(:,:,:,iFrame), [1 1], int2str(iFrame));
    end
end

%% Write video file
videoObject = VideoWriter(strcat(saveFolder, saveFile),'Uncompressed AVI');
videoObject.FrameRate = frameRate;
open(videoObject)
for iFrame = 1:length(frameNames)
    writeVideo(videoObject,videoFrames(:,:,:,iFrame))
end
close(videoObject)