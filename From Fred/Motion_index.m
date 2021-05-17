function [MI] = Motion_index(video, Timestamp, BG_frame)
% determine a motion index based on Jelitai et al., Nature communications,
% 2016.
% Note. HG 1/3/2017: Adapted from Fred's code.
%
% the video is a tif file. The time stamps of the video is used to
% determine the number of frame and is added to the matrix MI. Timestamp
% should be a one column matrix. BG_frame is the frame number used to substract the background.

% example: [MI] = Motion_index('example.tif', time, 100).

MI = zeros(size(Timestamp,1),1);

for i = 1:size(Timestamp,1)-1
    [X, ~] = imread(video, i);
    %X = x-BG;
    [X1, ~] = imread(video, i+1);
    %X1 = x1-BG;
    %temp = zeros(size(X,1), size(X,2));
%     for row = 1:size(X,1)
%         for col = 1:size(X,2)
%             temp(row,col) = X1(row,col) - X(row,col);
%         end
%     end
    FrameDiff = X1(:,:,1)-X(:,:,1);
    MI(i,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;
    %disp(i)
end

m = max(MI);
MI(:,1) = MI(:,1) ./ m; 

MI(:,2) = Timestamp(:,2);

end