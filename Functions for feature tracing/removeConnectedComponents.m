function [classification] = removeConnectedComponents...
    (predictorLocation, classification, sizeCut)
%Removes all positively identified pixels from classification that do
%not belong to connected components of at least the specified size as 
%evaluated by the locations in predcitorLocation.
%   Generates an image in which all positive classifications are = 1 and
%   all negative = 0. Then removes all positive classifications that belong
%   to connected components smaller than the specified sizeCut. This can
%   lead to no positive classifications remaining.
%
%   INPUT
%   predictorLocation are the y,x-coordinates of the upper left corner of
%       every section in the classification with the corresponding row.
%   classification denotes for the corresponding row in a predictor whether
%       it is a positive classification (1) or not (0).
%   sizeCut is the minimum size of the remaining connected components.
%
%   OUTPUT
%   classification denotes for the corresponding row in a predictor whether
%       it is a positive classification (1) or not (0). Only
%       classifications belonging to the largest component are still
%       classified as positive.

%% Validate input
validateattributes(predictorLocation,{'numeric'},{'positive', 'integer',...
    'size', [NaN,2]}, 'removeConnectedComponents', 'predictorLocation')
if iscell(classification)
    classification = cell2mat(classification);
end
validateattributes(classification,{'char'},{'ncols', 1}, ...
    'removeConnectedComponents', 'classification')
assert(length(classification) == size(predictorLocation,1), ...
    'predictorLocation and classification are of different length')
validateattributes(sizeCut,{'numeric'},{'positive', ...
    'integer', 'size', [1,1]}, 'removeConnectedComponents', 'sizeCut')

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
    numPixels = cellfun(@numel,connectedComponents.PixelIdxList);
    for iComponent = 1 : length(numPixels)
        if numPixels(iComponent) < sizeCut
            referenceFrame...
                (connectedComponents.PixelIdxList{iComponent}) = 0;
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