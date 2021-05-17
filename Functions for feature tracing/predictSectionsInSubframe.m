function [predictorLocation, classification, classificationScore] = ...
    predictSectionsInSubframe(model, imageFrame, sectionSize, ...
    subStartY, subStartX, subFrameSize)

%Computes a predictor from the subframe of the imageFrame using the 
%supplied model. The size of the subframe is specified by subStartY, 
%subStartX and subFrameSize. Then updates the predictorLocations to match 
%the original image.
%Important: sections of sectionSize will be predicted from the subframe,
%this should be considered when choosing the subFrameSize. Example: if the
%sectionSize is [25 25] and the subFrameSize [25 25] this will only result  
%in 1 prediction.
%
%INPUT
%   model = model generated from training data using fitcsvm, TreeBagger or
%       other compatible models. 
%   imageFrame = a two-dimensional array of non-negative numbers with
%       the image from which subframe should be extracted.
%   sectionSize = a 1x2 array of two non-zero natural numbers 
%      denoting the y,x size of the sections. sectionSize also needs to be 
%      smaller than the dimensions of the frame.
%   subStartY, subStartX = y- & x-coordinate respectively of the starting  
%       point (upper left corner) of the detail of the full frame for    
%       which the prediction is generated. If subStartY is smaller than 1 
%       it will automatically default to 1, so the OUTPUT parameters will
%       be shorter than expected.
%   subFrameSize = size of the subFrame from which the ouput parameters
%       will be predicted.
%
%OUTPUT
%   predictorLocation = the y,x-coordinates of the upper left corner of
%       every section in the predictor with the corresponding row.
%       Locations are matched to the original imageFrame.
%   classification = denotes for the corresponding row in the predictor
%       whether it is a positive classification (1) or not (0).
%   classificationScore = for each row in the predictor denotes the
%       probabilities to be classified into either 0 
%       (classificationScore(:,1)) or 1 (classificationScore(:,2)).

%% Validate INPUT (input model not validated)
validateattributes(imageFrame,{'numeric'},{'2d', 'integer', ...
    'nonnegative'}, 'predictSectionsInSubframe', 'imageFrame')
validateattributes(sectionSize,{'numeric'},{'positive', 'integer', ...
    'size', [1,2]}, 'predictSectionsInSubframe', 'sectionSize')
validateattributes(subStartY,{'numeric'},{'size', [1,1], 'integer', ...
    'positive'}, 'predictSectionsInSubframe', 'subStartY')
validateattributes(subStartX,{'numeric'},{'size', [1,1], 'integer', ...
    'positive'}, 'predictSectionsInSubframe', 'subStartX')
validateattributes(subFrameSize,{'numeric'},{'positive', 'integer', ...
    'size', [1,2]}, 'predictSectionsInSubframe', 'subFrameSize')

assert(min(size(imageFrame)-sectionSize) >= 0, 'sectionSize is too large')
assert(min(subFrameSize-sectionSize) >= 0, ...
    'sectionSize is too large or subFrameSize too small')
assert(min(size(imageFrame)-subFrameSize) >=0, 'subFrameSize is too large')


%% Assessing if the selected subframe is within the boundaries of the image
subEndY = subStartY + subFrameSize(:,1)-1;
subEndX = subStartX + subFrameSize(:,2)-1;

if subStartY < 1
    subStartY = 1;
end
if subStartX < 1
    subStartX = 1;
end
if subEndY > size(imageFrame,1)
    subEndY = size(imageFrame,1);
end
if subEndX > size(imageFrame,2)
    subEndX = size(imageFrame,2);
end

%% Generating subframe and doing prediction
subFrame = imageFrame(subStartY:subEndY, subStartX:subEndX);

[predictor, predictorLocation] = ...
    computePredictor(subFrame, sectionSize, 1);

[classification, classificationScore] = predict(model, predictor);

predictorLocation(:,1) = predictorLocation(:,1) + (subStartY-1);
predictorLocation(:,2) = predictorLocation(:,2) + (subStartX-1);
end