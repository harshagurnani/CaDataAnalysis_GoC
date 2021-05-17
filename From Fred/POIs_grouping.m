%% regroup POIs per axon.

%clear group Axon_dFF TimeAxon Numb_total_axon counter

group{1,1} = [1,2,3]; % the numbers are POI numbers.
 group{1,2} = [5];
 group{1,3} = [6,7,10,14,15,17,31,32,35];
 group{1,4} = [11];
 group{1,5} = [22,23];
 group{1,6} = [24,25,26,27,28,29];
 group{1,7} = [36,37,38,39,40,49,50,52,57,58,59,61,62,63,64,65,75,76,77,78,85,86,87];
 group{1,8} = [53,54,55,56];
 group{1,9} = [79];
 group{1,10} = [89,92,93,94,60];
 group{1,11} = [95,96,97];
 group{1,12} = [66,67,69];
% group{1,13} = [31,32,35]; 
% group{1,14} = [60]; 
% group{1,15} = [75,76,77,78,85,86,87];
% group{1,16} = [85,86,87];
% group{1,17} = [91];
% group{1,18} = [93];
% group{1,19} = [94];
% group{1,20} = [96];
% group{1,21} = [97];
% group{1,22} = [98];
% group{1,23} = [99:106,110:112];
% group{1,24} = [108];
% group{1,25} = [109];
% group{1,26} = [13:15];
%group{1,27} = [16:21];

for i = 1:size(group,2);
    Axon_dFF(:,i) = mean(dFF_POIs_data(:,group{1,i}),2);
    TimeAxon(:,i) = mean(POIs_data_time(:,group{1,i}),2);
end
Numb_total_axon = size(Axon_dFF,2);

s = 10; % smoothing factor.
[m, n] = size(Axon_dFF);
Data_smooth = smooth(Axon_dFF, s);
Axon_dFF_smooth = reshape(Data_smooth,m,n);

for i = 1:Numb_trials;
    if i == 1;
        Axon_dFF_smooth_trial{1,i} = Axon_dFF_smooth(1:size(dFF_POIs_data_trial{1,1},1),:);
        TimeAxon_trial{1,i} = TimeAxon(1:size(dFF_POIs_data_trial{1,1},1),:);
        counter = size(dFF_POIs_data_trial{1,1},1);
    end
    if i > 1;
        Axon_dFF_smooth_trial{1,i} = Axon_dFF_smooth(counter+1:(size(dFF_POIs_data_trial{1,1},1)*i),:);
        TimeAxon_trial{1,i} = TimeAxon(counter+1:(size(dFF_POIs_data_trial{1,1},1)*i),:);
        counter = size(dFF_POIs_data_trial{1,1},1)*i;
    end
end

% figure with traces for each axon.

y_offset = 1;
figure; hold on
for i = 1:size(Axon_dFF_smooth, 2);
    Y_offset = (y_offset * i) - y_offset;
    plot(TimeAxon(:,1)./1000, Axon_dFF_smooth(:,i) + Y_offset);
end
axis([0 inf -0.2 inf]);
xlabel('Time (s)');
ylabel('Number of axons');
title('dFF smooth all axons');
hold off
saveas(gcf, 'dFF smooth all axons.fig');

% figure tracing the axons.

figure; hold on
axon_vector_color = hsv(size(group,2));
POIs_numb = num2str((1:size(POI_coordinates.POIs_X,1))');
for i = 1:size(group,2);
    plot3(POI_coordinates.POIs_X(group{1,i},1), POI_coordinates.POIs_Y(group{1,i},1), POI_coordinates.POIs_Z(group{1,i},1), 'Color', axon_vector_color(i,:), 'LineWidth', 2);
end
view(3)
axis ij
hold off
hold on
for i = 1:size(group,2);
    plot3(POI_coordinates.POIs_X(group{1,i},1), POI_coordinates.POIs_Y(group{1,i},1), POI_coordinates.POIs_Z(group{1,i},1), 'o');
end
hold off
text(POI_coordinates.POIs_X, POI_coordinates.POIs_Y, POI_coordinates.POIs_Z, POIs_numb, 'FontSize', 10, 'HorizontalAlignment', 'right');
saveas(gcf, 'Axons tracing.fig');

