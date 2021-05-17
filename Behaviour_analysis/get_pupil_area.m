function [pupil] = get_pupil_area(video, threshDark, opts )
%% 



tic
%% If given vector

use_tiffile = true;
if use_tiffile, contents = dir(['*',video,'*.tif']); end
if size(contents)>1, warning('Too many pupil Tifs, using first one to setup threshold, rest will be processed iteratively'); end
if isempty(contents), use_tiffile = false; end
if ~use_tiffile,contents = dir(['*',video,'*.jpg']);
    n = natsortfiles({contents.name});
    [~,ndx] = natsortfiles({contents.name});
    contents = contents(ndx);
end


I = imread(contents(1).name);
picEyeAvg = mean(I,3);

% choose region of interest around eye
figure
uiwait(msgbox('Select region of interest. Double click to close loop. '));
maskEye = roipoly(picEyeAvg./255);
close all

%% Choose threshold for pupil value

DarkThreshFixed = false;
if isempty(threshDark), threshDark = 70; end

% check threshold value
figure
while ~DarkThreshFixed
uiwait(msgbox(sprintf('Choose threshDark value so that most pixels of pupil are showing. Current value is %d', threshDark) ));
for iTest = 1:10
    randFrame = ceil(rand*size(contents,1));
    disp(['random frame #' num2str(randFrame)])
    
    randImg = imread(contents(randFrame).name);
    randImg = randImg(:,:,1);
    
    randImg_threshed = randImg<threshDark;
    imshow(randImg_threshed.*maskEye)
    pause(1)
end
DarkThreshFixed = strcmpi(input('Was threshold value acceptable? Press y and Enter for Yes, or no to test again:','s'),'y');
if ~DarkThreshFixed, threshDark = input('Enter new threshold value:'); end
end
close all

pupil = nan(numel(n),1);
for k = 1:numel(n)
    filename(k) = {contents(k).name};
    I = imread(filename{k});
    if mod(k,500)==0
        disp([num2str(k) ' / ' num2str(size(n)) ' image ']);
    end
    
    % define image frame and threshold
    im_threshed = (I<threshDark) .* maskEye;
    im_threshed = imdilate( im_threshed, strel('square', opts(1)) );
%     im_threshed = imfill( im_threshed, 'holes' );
    pupil(k)  = sum(im_threshed(:));
%     [~, threshold] = edge(im, 'canny');
%     im_edges = edge(im, 'canny', threshold * opts(2)) .* maskEye;
    
end

end