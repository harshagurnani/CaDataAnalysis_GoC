%% import the timestamps for the eye- and whiskers-cam and few parameters on the video.

%% Whiskers video.

Whiskers.timestamp = importdata('WhiskersCam-relative times.txt'); % time stamp for the whisker camera.
Whiskers.TimeTotal = Whiskers.timestamp(end,2); % time in ms.
Whiskers.Numb_timestamp = size(Whiskers.timestamp,1);
Whiskers.Obj = VideoReader('WhiskersCam-1.avi');
% % Whiskers.Numb_video_frames = 0;
% % while hasFrame(Whiskers.Obj);
% %     readFrame(Whiskers.Obj);
% %     Whiskers.Numb_video_frames = Whiskers.Numb_video_frames + 1;
% % end
Whiskers.TimeTotal_video = Whiskers.Obj.Duration*1000; % time in ms.

% create a matrix with the video properties but too heavy for the whiskers
% cam.
% % Height = Obj.Height;
% % Width = Obj.Width;
% % Whiskers_video = zeros(Numb_video_frames, Height, Width);
% % for f = 1:Numb_video_frames;
% %     thisFrame = read(Obj, f);
% %     Whiskers_video(f,:,:) = rgb2gray(thisFrame);
% % end

%% Face cam (Eye cam) video.

FaceCam.timestamp = importdata('EyeCam-relative times.txt'); % time stamp for the face camera.
FaceCam.TimeTotal = FaceCam.timestamp(end,2); % time in ms.
FaceCam.Numb_timestamp = size(FaceCam.timestamp,1);
FaceCam.Obj = VideoReader('EyeCam-1.avi');
% % FaceCam.Numb_video_frames = 0;
% % while hasFrame(FaceCam.Obj);
% %     readFrame(FaceCam.Obj);
% %     FaceCam.Numb_video_frames = FaceCam.Numb_video_frames + 1;
% % end
FaceCam.TimeTotal_video = FaceCam.Obj.Duration*1000; % time in ms.


