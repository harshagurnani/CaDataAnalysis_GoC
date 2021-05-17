function [ contents, pupilFit ] = pupil_analysis( datapath, filename, code_loc, varargin )
%% To detect pupil: threshold, find boundary and fit ellipse. Returns structure with pupil's boundary points, centroid and axis lengths of fitted ellipse.
%
% A threshold is used to find all dark pixels, the largest connected dark
% component is kept (as the pupil), all holes are filled (for example to
% remove reflections inside the pupil, the edge of this object detected,
% and the edge points passed to fit_ellipse.m. The length of the short/long
% axis of the fitted ellipse, and the original centroid of the detected
% dark component are used as pupil size and location.
%
% Input:
% --------------------
%   - datapath      [STR]   Directory/path (absolute or relative) that
%                           contains image stack or sequence of eye cam
%                           videos. If empty, it searches in current
%                           directory.
%   - filename      [STR]   complete or partial filename - only load files
%                           that contain this fragment. Can be empty.
%   - code_loc      [STR]   Path to analysis codes. Can be empty if the
%                           analysis directory is already in path.
%   Optional:
%   - do_plot       [BOOL]  Whether to plot random images with pupil
%                           boundaries, and profile of long ans short axis
%                           changes in time. False by Default.
%                           Set to the optional arg that is of type 'logical' 
%
%   <    Threshold for Dark_pixels = thresh_factor * prctile( image, dark_prctile )    >
%   - thresh_factor [NUM]   Read as the first numeric optional argument
%   - dark_prctile  [NUM]   Read as the second (or more) numeric optional
%                           argument
%
% Output:
% ------------------
%   - files      [STRUCT ARRAY]     Contents of 'dir' - all files used in
%                                   analysis (helpful to know order of
%                                   frames if multiple files were used)
%   - pupilFit (STRUCT ARRAY) with fields:
%       - Boundary   [Cell array]   (x,y) coordinates of the points forming
%                                   the boundary of the detected pupil -
%                                   used as input to fit_ellipse.m
%       - CentroidX  [NUM]          x coordinate of detected pupil object
%       - CentroidY  [NUM]          y coordinate of detected pupil object
%       - Threshold  [Scalar unit8] Threshold value used to detect pixels
%                                   darker than this as corresponding to
%                                   the pupil
%       - LongAxis   [NUM]          Length of long axis at each frame/image
%       - ShortAxis  [NUM]          Length of short axis at each frame/image
%
%  Example Usage:
% -------------------------
%   [files, result] = pupil_analysis( '/Exp_1/', 'Mouse_HG001', 'D:/Silver-Lab-Project' ) 
%   <  Loads all files in '/Exp_1' containing 'Mouse_HG001' in its name.
%      Also allows scripts in D:/Silver-Lab-Project to be accessible. >
%
%   [files, result] = pupil_analysis( '/Exp_1/', 'Mouse_HG001', '\codes',  true, 1.2, 12 ) 
%   <  Loads all files in '/Exp_1' containing 'Mouse_HG001' in its name.
%      Also allows scripts in 'CurrentDirectory/codes' to be accessible.
%      Sets plotting to be true, threshold fudge factor to 1.2, and initial
%      threshold for dark pixels to be 12th percentile of pixel values in
%      image. >
%   Note that order of logical and numeric optional args is irrelevant i.e.
%   a similar result can be obtained as follows:
%   [files, result] = pupil_analysis( '/Exp_1/', 'Mouse_HG001', '\codes',  1.2, true, 12 ) %
%
%   However, exchanging the position of 1.2 and 12 will instead set
%   thresh_factor to 12, because it is always read as the first numeric
%   optional args, ni matter how many logical args come before it.
%
%
% PRE-REQUISITES:
% -------------------------
%       fit_ellipse.m   	 from https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
%                            Conic Ellipse representation = a*x^2+b*x*y+c*y^2+d*x+e*y+f=0
%                            (Tilt/orientation for the ellipse occurs when the term x*y exists (i.e. b ~= 0))
%
% HG 
% 30-1-2018

%% Parse arguments
if ~isempty(varargin)
    thresh_set = false;
    for ni = 1:numel(varargin)
        if islogical( varargin{ni} ),    do_plot = varargin{ni};    end     % Read logical argument as do_plot  
        if isnumeric( varargin{ni} )
           if thresh_set, dark_prctile  = varargin{ni};
           else,          thresh_factor = varargin{ni};     thresh_set = true;  end
        end
    end
end
% Set remaining to default
if ~exist( 'do_plot', 'var'),           do_plot       = false;  end
if ~exist( 'thresh_factor', 'var' ),    thresh_factor = 1.4;    end
if ~exist( 'dark_prctile', 'var'),      dark_prctile  = 10;     end

if isempty( datapath), datapath = ''; end
if isempty( code_loc), code_loc = ''; end

if isempty( filename), filename = ''; end   
filename = ['*', filename];    

%% Change to data directory/ load file in current directory 
% - all .jpg/.tif files matching filename or with specified extension

if isfolder(datapath),     cd location;  end
if contains( filename, '.' )
    contents = dir( filename );                                         % extension specified
else
    contents = [ dir([filename,'*.jpg']); dir([filename,'*.tif']) ];    % Load all jpg/tif files in datapath
end

% Directory where your analysis codes are :
if isfolder(code_loc), addpath( genpath(code_loc));   end

if isempty( contents)
    error( 'No matching files found')
end

if size(contents,1) == 1
    % single image stack 
    stack = true;
    allinfo = imfinfo(contents(1).name);        nFrames = numel(allinfo);
    contents = contents(1);
else
    % image sequence stored as separate files
    stack = false;
    nFrames = numel(contents);
end

%% Initialisation
[LongAxis, ShortAxis, EllipseAngle, A, B, X0, Y0] = deal(nan( nFrames, 1));
PupilCentroid   = nan(nFrames, 2);
PupilBoundary   = cell( nFrames,1);
Threshold       = nan(nFrames, 1);


test = imread( contents(1).name );  % read first frame
if size(test, 3) == 3, imgRGB = true; else, imgRGB = false; end

%% Run parallel on all frames
parfor frame = 1:nFrames
    if stack, curr_image = imread( contents.name, frame);
    else,     curr_image = imread( contents(frame).name); 
    end
    if imgRGB, curr_image = rgb2gray(curr_image);   end
    
    Threshold( frame )  = thresh_factor*  prctile( curr_image(:) , dark_prctile );
    pupil               = curr_image < Threshold(frame);
    pupil               = bwareafilt(pupil, 1);       % keep largest connected dark component
    pupil               = imfill( pupil, 'holes' );
    
    [ally, allx] = find( pupil );  
    PupilCentroid(frame,:) = [ mean( allx ), mean(ally)]; 
    
    tmp = edge(pupil);
    [ PupilBoundary{ frame }(:,2), PupilBoundary{ frame }(:,1)] =find( tmp ); 
    pupil_ellipse = fit_ellipse( PupilBoundary{frame}( :,1), PupilBoundary{frame}( :,2) );
    LongAxis( frame) = pupil_ellipse.long_axis;    ShortAxis(frame) = pupil_ellipse.short_axis;
%     A(frame)         = pupil_ellipse.a;            B(frame)         = pupil_ellipse.b;
%     X0(frame)        = pupil_ellipse.X0;           Y0(frame)        = pupil_ellipse.Y0;
    EllipseAngle(frame) = pupil_ellipse.phi;
end

pupilFit = struct(  'LongAxis',     mat2cell(LongAxis, nFrames), ...
                    'ShortAxis',    mat2cell(ShortAxis, nFrames), ...
                    'CentroidX',    mat2cell(PupilCentroid(:,1), nFrames), ...
                    'CentroidY',    mat2cell(PupilCentroid(:,2), nFrames), ...
                    'Boundary',     mat2cell(PupilBoundary, nFrames), ...
                    'Threshold',    mat2cell(Threshold, nFrames), ...
                    'EllipseAngle', mat2cell(EllipseAngle, nFrames)         );
%                     'X_A',          mat2cell(A, nFrames), ...
%                     'Y_B',          mat2cell(B, nFrames), ...
%                     'X0',           mat2cell(X0, nFrames), ...
%                     'Y0',           mat2cell(Y0, nFrames), ...



%% Plotting
if do_plot
    nPlot           = min( 16, nFrames );
    Frames_to_plot  = sort(randperm(nFrames, nPlot));
    
    figure;
    colormap(gray)
    nr  = ceil(sqrt(nPlot) );
    for seq = 1:nPlot
        subplot(nr,nr, seq)
        frame = Frames_to_plot(seq);
        if stack, curr_image = imread( contents.name, frame);
        else,     curr_image = imread( contents(frame).name);   end
        if imgRGB, curr_image = rgb2gray(curr_image);   end
        imagesc( curr_image ); hold on
        scatter( PupilBoundary{ frame }(:,1), PupilBoundary{ frame }(:,2), 6, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
        scatter( PupilCentroid( frame,1 ), PupilCentroid( frame,2 ), 'w', 'filled' )
%         ellipse( A(frame),B(frame),EllipseAngle(frame),X0(frame),Y0(frame),'b');
        title(['Frame ', num2str(frame)] )
    end
    
    
    figure;
    plot( LongAxis , 'r-' ); hold on
    plot( ShortAxis, 'b--' );
end