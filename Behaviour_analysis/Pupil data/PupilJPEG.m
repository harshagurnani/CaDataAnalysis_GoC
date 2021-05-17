%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README
%
% author: Jantine Broek
% e-mail: jantine.broek@yale.edu
% date: May 2017
% for: McCormick lab, Yale University, New Haven, USA
%
%
% Pre-req:
%            fit_ellipse  -  from https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
%                            Conic Ellipse representation = a*x^2+b*x*y+c*y^2+d*x+e*y+f=0
%                            (Tilt/orientation for the ellipse occurs when the term x*y exists (i.e. b ~= 0))
%                            EDIT: made changes in this function to create
%                            plots.
%
%            natsortfiles -  from http://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort
%
%            ellipse      -  from https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m.
%                            Only used to draw ellipse in figure.
%
%            inpaintn     - from https://www.mathworks.com/matlabcentral/fileexchange/27994-inpaint-over-missing-data-in-1-d--2-d--3-d--n-d-arrays
%                           To interpolate removed outliers
%
%
% Experimental set-up:
%            mice         - The mice are head-fixed and run on a wheel
%                           while recordings are taken.
%
%            camera       - Frame rate of the camera is 30Hz.
%
%            calibration  - An image of a ruler is taken with the camera in
%                           order to know the mm/pixel conversion. 
%
%            normalization- For every session, pupil measurements are normalized
%                           to the frame when the pupil is maximally
%                           dilated, which is always when the mouse is
%                           running. Pupil size = 1 correspond to its
%                           maximaly dilated state.
%
%            eye movement - This is minimal in head-fixed mice. 
%
% Goal:
%           (1) Using edge detection, the pupil is found and the 
%           coordinates of the pupil location will be used to estimate the 
%           function for (2) ellipse fitting. With the formula for ellipse 
%           fitting, the (3) pupil axis and location parameter are extracted
%           and will be used to compare the dynamics of the pupil fluctuation.
%
%
% Input:
%           .JPEG - these are files of individual pupil images obtained 
%                   with Virtual Dub software: 
%                   http://www.virtualdub.org/download.html.
%                   A sesion folder contains 10,000 frames, which is ~400MB
%                   Frame rate of the camera is 30Hz
%
%
% Manual input:
%           filepath - indicated in "Load Data" section
%           treshold - value for outlier detection, in "Outlier Removal"
%                      section
%
%
% Output:
%           longAxis   - size of the long axis of the ellipse
%           shortAxis  - size of the short axis of the ellipse
%           pupilXY    - pupil location
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
cd 'D:\myStuff\data_copy_from_server\HG 13 Emmy Noether\exp 02\180827_11_49_28 VidRec\split'
yourFolder = 'D:\myStuff\pupil data\Pupil data';
addpath(yourFolder);

% place files in natural order
contents = dir('*.jpg');
n = natsortfiles({contents.name});
[~,ndx] = natsortfiles({contents.name});
contents = contents(ndx);


%% Select eye ROI

% first frame
I = imread(contents(1).name);
picEyeAvg = mean(I,3);

% choose region of interest around eye
figure
uiwait(msgbox('Select region of interest. Double click to close loop. '));
maskEye = roipoly(picEyeAvg./255);
close all

%% Choose threshold for pupil value

% threshold (default is 25)
threshDark = 11;

% check threshold value
figure
uiwait(msgbox('Choose threshDark value so that most pixels of pupil are showing.'));
for iTest = 1:15
    randFrame = ceil(rand*size(contents,1));
    disp(['random frame #' num2str(randFrame)])
    
    randImg = imread(contents(randFrame).name);
    randImg = randImg(:,:,1);
    
    randImg_threshed = randImg<threshDark;
    imshow(randImg_threshed.*maskEye)
    pause
end
close all


%% Edge detection and Ellipse fitting

% eye may be closed if very few dark pixels found (default = 200)
threshNdarkPix = 200;

% fudge factor for estimated radius (default = 1.4)
radiusFudge = 1.4;

% canny edge detection threshold fudge factor (default 1). Lower means more edges will be found
edgeFudge = 1;

% edge pixels must be at least this close to dark pixels (default 3)
threshPupilDist = 3;

% edge pixels must be continguous with at least this many other pixels (default 15)
minEdgeSize = 10;

% preallocate variables
filename = cell(1,numel(n));
shortAxis = zeros(1, numel(n));
longAxis = zeros(1, numel(n));
pupilXY = zeros(1, numel(n));

tic
for k = 1:numel(n)
    filename(k) = {contents(k).name};
    I = imread(filename{k});
    
    % create title for images
    [pathstr, name, ext] = fileparts(filename{k});
    ix=[ numel(name)+1, strfind(name,'_')];     %Name may not have underscores or 4 underscores!!!!
    t = name(1:ix(end)-1);                      
    
    
    % display loop progression
    if mod(k,500)==0
        disp([num2str(k) ' / ' num2str(size(n)) ' image ']);
    end
    
    
    % define image frame and threshold
    im_threshed = (I<threshDark) .* maskEye;
    
    % estimate pupil center and size based of number of pixels
    [row1, col1] = find(im_threshed);
    
    
    % check to see if eye is open
    if length(row1) > threshNdarkPix
        
        % pupil center
        pupilCenter_estimate = [median(row1), median(col1)];
        
        % pupil radius (from A=pi*r^2)
        pupilRadius_estimate = sqrt(length(row1)/pi);
        
        
        % delete pixels which are too far away from pupil
        for l = 1:length(row1)
            
            % distance from estimated center
            dXY = pupilCenter_estimate - [row1(l), col1(l)];
            d = sqrt(sum(dXY .* dXY));
            
            % delete pixels if too far
            if d > pupilRadius_estimate * radiusFudge
                im_threshed(row1(l), col1(l)) = false;
            end
            
        end
        
        
        % canny edge detection
        if size(I,3)==3, im = rgb2gray(I); else, im = I; end        %HG: Note: Image may already be grayscale
        [~, threshold] = edge(im, 'canny');
        im_edges = edge(im, 'canny', threshold * edgeFudge) .* maskEye;
        
        % delete edge pixels which are too far from center
        [row2, col2] = find(im_edges);
        for m = 1:length(row2)
            
            % distance from estimated center
            dXY = pupilCenter_estimate - [row2(m), col2(m)];
            d = sqrt(sum(dXY .* dXY));
            
            % delete edge pixels which are too far from dark pixels
            [~,d] = dsearchn([row1, col1], [row2(m), col2(m)]);
            if d > threshPupilDist
                im_edges(row2(m),col2(m)) = false;
            end
            
        end
        
        
        % delete edge pixels which are not part of large group
        % label ROIs by unique number
        im_edges2 = bwlabel(im_edges);
        
        % number of ROIs in edges matrix
        nRois = max(im_edges2(:));
        
        % clear edge matrix and add back in if ROI is big
        im_edges = false(size(im_edges));
        for iROI = 1:nRois
            roiSize = sum(sum(im_edges2 == iROI));
            
            if roiSize >= minEdgeSize
                im_edges(im_edges2==iROI) = true;
            end
            
        end
        
        % ellipse fitting if there are sufficient edges pixels
        if sum(im_edges(:)) > minEdgeSize
            
            % fit ellipse
            [row3, col3] = find(im_edges);
            pupilEllipse = fit_ellipse(row3, col3);
            
            % if ellipse found, use long axis as pupil diameter
            if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
                
                % record pupil axis
                shortAxis(k) = pupilEllipse.short_axis;
                longAxis(k) = pupilEllipse.long_axis;
                
                % record pupil position
                pupilXY(1,k) = pupilEllipse.Y0_in;
                pupilXY(2,k) = pupilEllipse.X0_in;
                
            end
            
        end
        
    end
    
    
    
    %   % plot pic with fitted ellipse
    %     figure;
    %     imshow(im); title(sprintf('Pupil Ellipse 2, image %s', t));
    %
    %     % resize figure
    %     set(gcf, 'Position', [150, 150, 400, 400]);
    %
    %     if ~isempty(pupilEllipse)
    %
    %         ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
    %
    %         % pupil center
    %         hold on
    %         scatter(pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
    %
    %         % pupil axis
    %         cos_phi = pupilEllipse.cos_phi;
    %         sin_phi = pupilEllipse.sin_phi;
    %
    %         % rotation matrix to rotate the axes with respect to an angle phi
    %         R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    %
    %         % the axes
    %         ver_line        = [ [pupilEllipse.X0 pupilEllipse.X0]; pupilEllipse.Y0+pupilEllipse.b*[-1 1] ];
    %         horz_line       = [ pupilEllipse.X0+pupilEllipse.a*[-1 1]; [pupilEllipse.Y0 pupilEllipse.Y0] ];
    %         new_ver_line    = R * ver_line;
    %         new_horz_line   = R * horz_line;
    %
    %         % draw major axis
    %         plot( new_ver_line(2,:),new_ver_line(1,:),'b' );
    %         % draw minor axis
    %         plot( new_horz_line(2,:),new_horz_line(1,:),'g' );
    %
    %         hold off
    %      end
    
end
toc

%% Extract tuple variables for DataJoint
% extract session date, session number, filename for x-axis, movie number (first number of filenr), frame
% number (rest of numbers in filenr), animal id
token = strtok(filename,'.');
D = regexp(token, '_', 'split');
D = vertcat(D{:});

% filenr  - frame number (either the last number or the 5th one
if size(D,2)>4, filenr_temp = D(:,5); else, filenr_temp = D(:,end);  end
G = sprintf('%s*', filenr_temp{:}); % change cell to double
filenr = sscanf(G, '%f*');
filenr = filenr';

% session date: should be in string format "2017-05-12"
if size(D,2)>4
    session_date_temp = D(:,1); 
    session_date = datetime(session_date_temp, 'InputFormat', 'ddMMyy', 'Format', 'yyyy-MM-dd');
%session_date = datetime(session_date_temp, 'InputFormat', 'MMddyy', 'Format', 'yyyy-MM-dd');
    session_date = string(session_date); % as string
else
    
end

% session number
if size(D,2)>4
    session_number_temp = D(:,3);
    F = sprintf('%s*', session_number_temp{:}); % change cell to double
    session_number = sscanf(F, '%f*');
else
    
    
end


% movie number and frame number
movie_number = zeros(1,numel(filenr));
frame_number = zeros(1,numel(filenr));
for m = 1:numel(filenr)
    % file number is currently only frame number, and no movie number. so
    % putting movie number as 1.
%     G = num2str(filenr(m));
%     movie_number(m) = str2double(G(1));
%     frame_number(m) = str2double(G(2:end));
    G = num2str(filenr(m));
    movie_number(m) = 1;%str2double(G(1));
    frame_number(m) = str2double(G);
end

% extract animal_id column
% animal_id_temp = D(:,2);
% animal_id_temp = animal_id_temp(:,1);
% animal_id_temp = sprintf('%s*', animal_id_temp{:});
% animal_id = regexp(animal_id_temp, '\d*', 'Match');
animal_id = 'FRED_TEST_ANIMAL';%str2double(animal_id');


%% Plot axis

figure;
hold on

subplot(2,1,1);
plot(shortAxis, '-.ob')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('shortAxis');

subplot(2,1,2);
plot(longAxis, '-.or')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('longAxis');
hold off


%% Collect raw data

data = [filenr; shortAxis; longAxis; pupilXY];
data(data == 0) = NaN;
data = data';

% find NaN
[row, col] = find(isnan(data));
nrNaN = unique(row);
dataNaN = data(nrNaN,:);
% check NaN images (optional to view them)
filename_NaN = filename(:,nrNaN)';
% figure;
% for ii = 1:numel(filename_NaN)
%    INaN = imread(filename_NaN{ii});
%    imshow(filename_NaN{ii})
%    pause;
% end


%% Outlier removal

% change 0 values to NaN
shortAxis(shortAxis == 0) = nan;
longAxis(longAxis == 0) = nan;
pupilXY(pupilXY == 0) = nan;

% compute the mean value
shortAxis_mean = nanmean(shortAxis);
longAxis_mean = nanmean(longAxis);
pupilXY_mean = nanmean(pupilXY,2);

% compute the absolute difference
shortAxis_absdiff = abs(shortAxis - shortAxis_mean);
longAxis_absdiff = abs(longAxis - longAxis_mean);
pupilXY_absdiff = abs(pupilXY - pupilXY_mean);

% compute the median of the absolute difference
shortAxis_mad = nanmedian(shortAxis_absdiff);
longAxis_mad = nanmedian(longAxis_absdiff);
pupilX_mad = nanmedian(pupilXY_absdiff(1,:));
pupilY_mad = nanmedian(pupilXY_absdiff(2,:));

% outliers if the absolute difference is moe than some factor times the mad value
% shortAxis
sensitivityFactor = 5; % change this to a lower value for more stringent outlier removal (default 5)
thresholdValue_sA = sensitivityFactor * shortAxis_mad;
outlierIndexes_sA = abs(shortAxis_absdiff) > thresholdValue_sA;
% longAxis
thresholdValue_lA = sensitivityFactor * longAxis_mad;
outlierIndexes_lA = abs(longAxis_absdiff) > thresholdValue_lA;
% pupilXY
thresholdValue_X = sensitivityFactor * pupilX_mad;
thresholdValue_Y = sensitivityFactor * pupilY_mad;
outlierIndexes_X = abs(pupilXY_absdiff(1,:)) > thresholdValue_X;
outlierIndexes_Y = abs(pupilXY_absdiff(2,:)) > thresholdValue_Y;

% extract outlier values
shortAxis_outliers = shortAxis(outlierIndexes_sA);
longAxis_outliers = longAxis(outlierIndexes_lA);
pupilX = pupilXY(1,:);
pupilY = pupilXY(2,:);
pupilX_outliers = pupilX(outlierIndexes_X);
pupilY_outliers = pupilY(outlierIndexes_Y);

% remove outliers from original dataset and replace with NaN
shortAxis(outlierIndexes_sA) = nan;
longAxis(outlierIndexes_lA) = nan;
pupilX(outlierIndexes_X) = nan;
pupilY(outlierIndexes_Y) = nan;


%% Interpolate NaN values

shortAxis_inter = inpaintn(shortAxis);
longAxis_inter = inpaintn(longAxis);
pupilX_inter = inpaintn(pupilX);
pupilY_inter = inpaintn(pupilY);

pupilXY_inter = zeros(2, numel(n)); 
pupilXY_inter(1,:) = pupilX_inter;
pupilXY_inter(2,:) = pupilY_inter;

% plot axis with interpolated values
figure;
plot(data(:,3), '-or')
hold on
plot(longAxis_inter, '-og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
title('Original and interpolated pupil long axis')
legend('Original', 'Outlier removal');
hold off

% plot pupilXY with interpolated values
figure;
hold on
subplot(2,1,1);
plot(data(:,4), '-or')
hold on
plot(pupilX_inter, '-og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('pupil location (px)');
title('Original and interpolated pupil location')
legend('Original', 'Outlier removal');

subplot(2,1,2);
plot(data(:,5), '-or')
hold on
plot(pupilY_inter, '-og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('pupil location (px)');
title('Original and interpolated pupil location')
legend('Original', 'Outlier removal');
hold off


%% Collect data

% store data with headers
coldata = {'Filenr', 'shortAxis', 'longAxis', 'Pupil_x', 'Pupil_y'};
data = array2table(data, 'VariableNames', coldata);

data_inter = [filenr; shortAxis_inter; longAxis_inter; pupilXY_inter]';
data_inter = array2table(data_inter, 'VariableNames', coldata);

% save raw data
fname = sprintf('dataRaw_Mouse%d.mat', animal_id(1:1));
save(fname, 'data');

% save interpolated data
fname = sprintf('dataInter_Mouse%d.mat', animal_id(1:1));
save(fname, 'data_inter');


%% Add to DataJoint

% addpath /media/jantine/Data/04_DataJoint/2PE/schemas
% 
% 
% animal_id = num2cell(animal_id);
% %session_date = num2cell(session_date);
% session_date = arrayfun(@char, session_date, 'uni', false);
% session_number = num2cell(session_number);
% movie_number = num2cell(movie_number');
% frame_number = num2cell(frame_number');
% short_axis = num2cell(shortAxis_inter');
% long_axis = num2cell(longAxis_inter');
% pupilXY_inter = pupilXY_inter';
% pupil_x = num2cell(pupilXY_inter(:,1));
% pupil_y = num2cell(pupilXY_inter(:,2));
% 
% tuple = horzcat(animal_id, session_date, session_number, movie_number, frame_number, short_axis, long_axis, pupil_x, pupil_y);
% 
% 
% % order structure as DataJoint keys
% a = pupil.EyeROI;
% fields = a.header.names;
% fields = fields(1:end-1);
% tuple = cell2struct(tuple, fields, 2);
% 
% 
% a.insert(tuple)













