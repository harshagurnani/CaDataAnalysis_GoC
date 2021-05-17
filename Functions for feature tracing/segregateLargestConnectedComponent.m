function [classification] = segregateLargestConnectedComponent...
    (predictorLocation, classification, varargin)
%Removes all positively identified pixels from classification that do
%not belong to the largest connected component (or the nth largest
%connected component specified) as evaluated by the locations in
%predcitorLocation.
%   Generates an image in which all positive classifications are = 1 and
%   all negative = 0. Unless remainingComponents is specified it then finds
%   the connected components and removes all but the biggest one. If
%   remainingComponents is specified it removes all but the specified
%   number of largest components. Example: remainingComponents = 3 will
%   result in the 3 largest components remaining if possible.
%
%   INPUT
%   predictorLocation are the y,x-coordinates of the upper left corner of
%       every section in the classification with the corresponding row.
%   classification denotes for the corresponding row in a predictor whether
%       it is a positive classification (1) or not (0).
%
%   OPTIONAL INPUT - pass as a structure with the appropriate fields
%   remainingComponents is the number of largest components remaining after
%       the removal of smaller components. Default = 1.
%
%   OUTPUT
%   classification denotes for the corresponding row in a predictor whether
%       it is a positive classification (1) or not (0). Only
%       classifications belonging to the largest component are still
%       classified as positive.

%% Validate input
validateattributes(predictorLocation, {'numeric'}, {'positive', 'integer',...
    'size', [NaN,2]}, 'segregateLargestConnectedComponent', ...
    'predictorLocation')
if iscell(classification)
    classification = cell2mat(classification);
end
validateattributes(classification,{'char'},{'ncols', 1}, ...
    'segregateLargestConnectedComponent', 'classification')
assert(length(classification) == size(predictorLocation,1), ...
    'predictorLocation and classification are of different length')
remainingComponents = 1;
if ~isempty(varargin)
    optionalInput = varargin{1};
    if isfield(optionalInput, 'remainingComponents')
        remainingComponents = optionalInput.remainingComponents;
        validateattributes(remainingComponents,{'numeric'},{'positive', ...
            'integer', 'size', [1,1]}, ...
            'segregateLargestConnectedComponent', 'remainingComponents')
    end
end

if ~isempty(classification == '1')
    %% Generate reference frame
    referenceFrame = zeros(max(predictorLocation(:, 1)), ...
        max(predictorLocation(:, 2)), 'uint8');
    for iLocation = 1:length(classification)
        if classification(iLocation) == '1'
            referenceFrame(predictorLocation(iLocation,1), ...
                predictorLocation(iLocation, 2)) = 1;
        end
    end
    %% Compute connected components
    connectedComponents = bwconncomp(referenceFrame,4);
    if connectedComponents.NumObjects > remainingComponents
        numPixels = cellfun(@numel,connectedComponents.PixelIdxList);
        for iComponent = 1 : length(numPixels)
            sortedNumPixels = sort(numPixels);
            if numPixels(iComponent) < ...
                    sortedNumPixels(end-remainingComponents+1)
                referenceFrame...
                    (connectedComponents.PixelIdxList{iComponent}) = 0;
            end
        end
    end
    %% Remove smaller connectedComponents from classification
    for iLocation = 1:length(classification)
        if classification(iLocation) == '1'
            referencePixel = referenceFrame...
                (predictorLocation(iLocation,1), ...
                predictorLocation(iLocation, 2));
            if referencePixel == 0
                classification(iLocation) = '0';
            end
        end
    end
end