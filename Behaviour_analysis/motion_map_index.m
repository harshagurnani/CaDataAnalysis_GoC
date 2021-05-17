function [MMap, MI] = motion_map_index( video, times, alpha, beta )
%% Based on Powell, Duguid et al. eLife (2015)
% Pubmed: <http://www.ncbi.nlm.nih.gov/pubmed/26083712 Pubmed link> 
% Motion map AND motion index
% SYNTAX: [mmap, mi] = motion_map_index2( 'C:/Harsha/Npad.tif', Time1, 0.3, 0.9 )
%%%%%%%%%%%%%
%%%% NOTE THAT MOTION MAP MAY BE OF A LARGE SIZE - better to do on
%%%% smaller ROI.
%%%%%%%%%%%%%
%
%
% Steps:
% Each time/frame has its own background(BG), which is a weighted average
% of the previous timestep's background and the current frame.
%   BG[x,y, t] = alpha * Frame[x,y,t] + (1 - alpha) * BG[x,y,t-1]
% Then an absolute difference image is computed between frame and
% background:
%   DiffIm[x,y,t] = | Frame[x,y,t] - BG[x,y,t] |
% Effectively this is an average of frame and previous background.
%
% Then, they get Motion Map by smoothening - computing weighted average of
% successive frames.
%   MMap[x,y,t] = beta * DiffIm[x,y,t] + (1-beta) * MMap[x,y,t-1]
%
% Motion index is just the number of pixels in the frame with motion map
% MMap > threshold (40 by default).
%
% Harsha Gurnani. March 2017

 %% No. of frames/Frame size
 nt = size(times, 1)-1;
 BG = im2double(imread(video,1));   %Background
 nx = size(BG,1);
 ny = size(BG,2);
 
 MMap = nan(nx,ny,2);
 MI = nan(nt,1);
 
 BG = alpha * im2double(imread(video, 2)) + (1-alpha) * BG;  %Updated background
 MMap(:,:,1) = abs( im2double(imread(video,2)) - BG );       %MMap = DiffIm for first step
 MI(1) = sum(sum(MMap(:,:,1)>=40/255 ));
 
 for t = 2:nt
     %Update background
     BG = alpha * im2double(imread(video, t+1)) + (1-alpha) * BG;  
     % |Frame - Background|
     DiffIm = reshape(abs( im2double(imread(video,t+1)) - BG ), [nx,ny,1]); 
     % Weighted average of Motion Map
     MMap(:,:,t) = beta * DiffIm + (1-beta) * MMap(:,:,t-1);
     % Get Motion Index: How many pixels had "motion"?
     MI(t) = sum(sum(MMap(:,:,2)>=40/255 ));
 end

 MaxM = prctile(MI,99.95);
 MinM = min(MI);%Minimum will always be zero.
 MI = MI/MaxM;
end