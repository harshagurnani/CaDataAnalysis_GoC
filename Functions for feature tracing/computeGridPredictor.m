function [predictor] = computeGridPredictor(gridSpacing, gridSize, ...
    centerOrigin, imageFrame, sectionSize)
%Returns the predictor of a grid of sections sourrounding the section with
%the centerOrigin.
%   Computes a grid of gridSize with a pixel spacing of gridSpacing around
%   the centerOrigin supplied. Then passes all these sections into a
%   predictor array where each row corresponds to one section. Sections
%   that are on the grid but not fully within the imageFrame are excluded,
%   therefore the length of the predictor might be shorter than gridSize^2.
%
%   INPUT:
%   gridSpacing must be a non-zero natural number larger than 1 and is the
%       pixel distance between the section in the grid of negative
%       sections.
%   gridSize must be a non-zero natural number larger than 1. It is the
%       one dimensional size of the grid. The number of sections on the
%       grid are gridSize^2.
%   centerOrigin must be a 1x2 array of two non-zero natural numbers
%      denoting the start (upper left corner) of the center section.
%   imageFrame must be a two-dimensional array of non-negative numbers with
%       the image from which sections should be extracted.
%   sectionSize must be a 1x2 array of two non-zero natural numbers
%      denoting the y,x size of the sections. sectionSize also needs to be
%      smaller than the dimensions of the frame.
%
%   OUTPUT:
%   predictor is a matrix in which each row (y dimension) is a linearized
%       section belonging to the grid.

%% Validate input
validateattributes(gridSpacing, {'numeric'}, ...
    {'positive','integer','>', 1}, 'computeGridPredictor', 'gridSpacing')
validateattributes(gridSize, {'numeric'}, {'positive', 'integer'}, ...
    'computeGridPredictor', 'gridSize')
validateattributes(centerOrigin,{'numeric'}, ...
    {'integer', 'size', [1,2]}, 'computeGridPredictor', 'centerOrigin')
validateattributes(imageFrame, {'numeric'}, ...
    {'2d', 'integer', 'nonnegative'}, 'computeGridPredictor', 'imageFrame')
validateattributes(sectionSize,{'numeric'},{'positive', 'integer', ...
    'size', [1,2]}, 'computeGridPredictor', 'sectionSize')
assert(min(size(imageFrame)-sectionSize) >= 0, 'sectionSize is too large')

%% Compute the grid coordinates
maxGridDistance = (gridSize-1)*gridSpacing;
yCoordinates = 0 : gridSpacing : maxGridDistance;
yCoordinates = yCoordinates + centerOrigin(1) - floor(maxGridDistance/2);
yCoordinates = repmat(yCoordinates,1,gridSize);
xCoordinates = 0 : gridSpacing : maxGridDistance;
xCoordinates = xCoordinates + centerOrigin(2) - floor(maxGridDistance/2);
xCoordinates = sort(repmat(xCoordinates,1,gridSize));
gridOrigins = vertcat(yCoordinates, xCoordinates);
if ismember(centerOrigin, transpose(gridOrigins), 'rows')
    [~, centerInGrid] = ismember(centerOrigin, ...
        transpose(gridOrigins), 'rows');
    gridOrigins(:, centerInGrid) = [];
end

%% Remove coordinates belonging to sections not fully inside the image
gridOrigins(:, gridOrigins(1,:) < 1) = [];
gridOrigins(:, gridOrigins(2,:) < 1) = [];
gridEndPoints = gridOrigins + ...
    transpose(repmat(sectionSize, size(gridOrigins,2), 1)) - 1;
imageSize = size(imageFrame);
gridOrigins(:, gridEndPoints(1,:) > imageSize(1,1)) = [];
gridEndPoints(:, gridEndPoints(1,:) > imageSize(1,1)) = [];
gridOrigins(:, gridEndPoints(2,:) > imageSize(1,2)) = [];
gridEndPoints(:, gridEndPoints(2,:) > imageSize(1,2)) = [];

%% Compute predictor of sections belonging to gridCoordinates
if ~isempty(gridOrigins)
    predictor = zeros(size(gridOrigins,2), prod(sectionSize), ...
        'like', imageFrame);
    for iGrid = 1:size(gridOrigins,2)
        section =imageFrame(gridOrigins(1,iGrid):gridEndPoints(1,iGrid),...
            gridOrigins(2,iGrid):gridEndPoints(2,iGrid));
        predictor(iGrid, :) = section(:);
    end
else
    predictor = zeros(1, prod(sectionSize), 'like', imageFrame);
end
end