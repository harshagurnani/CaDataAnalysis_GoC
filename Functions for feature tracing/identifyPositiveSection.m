function [section, sectionOrigin] = identifyPositiveSection(imageFrame, ...
    boundingBox, sectionSize)
% Returns the predictor of the image section in the imageFrame defined by
% the boundingBox and the sectionSize.
%
%   INPUT
%   imageFrame must be a two-dimensional array of non-negative numbers with
%       the image from which sections should be extracted.
%   boudingBox must be a 1x4 array of non-zero natural numbers denoting the
%      [x start, y start, width, heigth] of the boundingBox.
%   sectionSize must be a 1x2 array of two non-zero natural numbers
%      denoting the y,x size of the sections. sectionSize also needs to be
%      smaller than the dimensions of the imageFrame.
%
%   OUTPUT
%   section is the 2D matrix of the section of the imageFrame specified by
%       boundingBox und sectionSize.
%   sectionOrigin is the y,x-coordinate of the upper left corner of the
%       section.

%% Validate inputs
validateattributes(imageFrame,{'numeric'},{'2d', 'integer', ...
    'nonnegative'}, 'identifyPositiveSection', 'imageFrame')
validateattributes(boundingBox,{'numeric'},{'positive', 'integer', ...
    'size', [1,4]}, 'identifyPositiveSection', 'boundingBox')
validateattributes(sectionSize,{'numeric'},{'positive', 'integer', ...
    'size', [1,2]}, 'identifyPositiveSection', 'sectionSize')
assert(min(size(imageFrame)-sectionSize) >= 0, 'sectionSize is too large')

%% Identify Section

heightCorrection = floor((sectionSize(1,1)-boundingBox(4))/2);
widthCorrection = floor((sectionSize(1,2)-boundingBox(3))/2);
originY = boundingBox(2)-heightCorrection;
originX = boundingBox(1)-widthCorrection;

if min(originY, originX)>0 && ...
        min([originY originX]+sectionSize <= size(imageFrame))
    section = imageFrame(originY:originY+sectionSize(1,1)-1, ...
        originX:originX+sectionSize(1,2)-1);
else
    section = zeros(sectionSize);
end

sectionOrigin = [originY originX];

end