function [predictor, predictorLocation] = computePredictor(imageFrame, ...
    sectionSize, varargin)
% Computes a predictor matrix of sections of a supplied image frame.
% Returns the y,x origin (start) of every section in predictor.
%
%   The computed predictor matrix stores one image section in each row. If
%   sectionSpacing is set to one the shift between the sections will be one
%   pixel, if sectionSpacing is larger than one the shift will be the
%   corresponding number of pixels.
%
%   INPUT:
%   imageFrame must be a two-dimensional array of non-negative numbers with
%       the image from which sections should be extracted.
%   sectionSize must be a 1x2 array of two non-zero natural numbers
%      denoting the y,x size of the sections. sectionSize also needs to be
%      smaller than the dimensions of the frame.
%
%   OPTIONAL INPUT:
%   sectionSpacing must be a non-zero natural number. It is the pixel
%       distance in y or x between the sections. sectionSpacing = 1 will
%       return every possible section of the specified size in the image.
%       Default = 1.
%
%   OUTPUT:
%   predictor is a matrix in which each row (y dimension) is a linearized
%       section of the imageFrame.
%   predictorLocation are the y,x-coordinates of the upper left corner of
%       every section in the predictor with the corresponding row.

%% Validate the inputs
validateattributes(imageFrame, {'numeric'}, ...
    {'2d', 'integer', 'nonnegative'}, 'computePredictor', 'imageFrame')
validateattributes(sectionSize, {'numeric'}, {'positive', 'integer', ...
    'size', [1,2]}, 'computePredictor', 'sectionSize')
assert(min(size(imageFrame)-sectionSize) >= 0, 'sectionSize is too large')

if ~isempty(varargin)
    sectionSpacing = varargin{1};
    validateattributes(sectionSpacing, {'numeric'}, ...
        {'positive', 'integer'}, 'computePredictor', 'sectionSpacing')
else
    sectionSpacing = 1;
end

%% Compute the predictor
[imageSizeY, imageSizeX] = size(imageFrame);
nShiftsY = floor((imageSizeY-sectionSize(1,1))/sectionSpacing);
nShiftsX = floor((imageSizeX-sectionSize(1,2))/sectionSpacing);

predictor = zeros(nShiftsY*nShiftsX,prod(sectionSize), 'like', imageFrame);
predictorLocation = zeros(nShiftsY*nShiftsX, 2);

for iShiftY = 0:nShiftsY
    for iShiftX = 0:nShiftsX
        yOrigin = 1 + iShiftY*sectionSpacing;
        xOrigin = 1 + iShiftX*sectionSpacing;
        section = imageFrame(yOrigin:yOrigin+sectionSize(1,1)-1, ...
            xOrigin:xOrigin+sectionSize(1,2)-1);
        sectionCount = iShiftY*(nShiftsX+1) + iShiftX+1;
        predictor(sectionCount, :) = section(:);
        predictorLocation(sectionCount,1) = yOrigin;
        predictorLocation(sectionCount,2) = xOrigin;
    end
end

predictor = single(predictor);

end

