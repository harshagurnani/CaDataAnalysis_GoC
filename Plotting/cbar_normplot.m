function [] = plot_scaled_image( data, tp_roi_trial, ftype)
%% Given: 'data' with dim index array 'tp_roi_trial' ([tp, roi, trial]) and scaling type, plot scaled vetdion using imagesc with separate rows for each POI and trials one after the other along the row
temp = permute(norm,[1,3,2]);
Newtemp=[];
trial = size(norm,2);
poi=size(norm,3);
for tn = 1:trial
    Newtemp = vertcat(Newtemp,temp(:,:,tn));
    %Newtemp = [Newtemp; zeros(100,poi)];
end
clear temp
Newtemp = transpose(Newtemp);
clims = [0 1];
imagesc(Newtemp, clims)
