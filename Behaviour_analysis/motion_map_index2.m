function [MI] = motion_map_index2( video, times, alpha, beta, thresh )
%% Based on Powell, Duguid et al. eLife (2015)
% Pubmed: <http://www.ncbi.nlm.nih.gov/pubmed/26083712 Pubmed link> 
% Motion index based on motion map. 
%%%%%%%%%% This function does NOT return motion map.
% SYNTAX: 
%   mi = motion_map_index2( 'C:/Harsha/Npad.tif', timestamps, 0.3, 0.9, 30 )
%
% Analysis steps: 
% 1) Each time/frame has its own background(BG), which is a weighted average
% of the previous timestep's background and the current frame.
%   BG[x,y, t] = alpha * Frame[x,y,t] + (1 - alpha) * BG[x,y,t-1]
%
% 2) Then an absolute difference image is computed between frame and
% background:
%   DiffIm[x,y,t] = | Frame[x,y,t] - BG[x,y,t] |
%
% 3) Then, they get Motion Map by computing weighted average of
% successive frames' motion.
%   MMap[x,y,t] = beta * DiffIm[x,y,t] + (1-beta) * MMap[x,y,t-1]
%
% 4) Motion index is just the number of pixels in the frame with motion map
% MMap > threshold (40 by default).
%
% Harsha Gurnani. March 2017
tic
 %% No. of frames/Frame size
 nt = size(times, 1)-1;
 BG = im2double(imread(video,1));
 nx = size(BG,1);
 ny = size(BG,2);
 
 if ~(thresh >=1 && thresh<=255)
     thresh = 40;   %Threshold for change in pixel intensity (between 1 to 255)
     %Default = 40 as used in Powell, Duguid et al 2015
 end 
 
 MMap_curr = nan(nx,ny);
 MI = nan(nt,1);
 
 %Initial frame
 BG = alpha *im2double(imread(video, 2)) + (1-alpha) * BG;         %Updated background
 MMap_curr = abs( im2double(imread(video,2)) - BG );    %MMap = DiffIm for first step
 MI(1) = sum(sum( MMap_curr >= thresh/255 ));
 MMap_prev = MMap_curr;
 
 for t = 2:nt
     %Update background
     BG = alpha * im2double(imread(video, t+1)) + (1-alpha) * BG;      
     % |Frame - Background|
     DiffIm = abs( im2double(imread(video,t+1)) - BG );    
     % Weighted average of Motion Map
     MMap_curr = beta * DiffIm + (1-beta) * MMap_prev;
     MMap_prev = MMap_curr;
     % Get Motion Index: How many pixels had "motion"?
     MI(t) = sum(sum(MMap_curr >= thresh/255 ));
 end

 MaxM = prctile(MI,99.95);
 %Minimum will always be zero.
 MI = MI/MaxM;
 toc
end