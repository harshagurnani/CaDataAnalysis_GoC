% updated 02/12/2016 with time between trial incorporated in the POI time
% stamps. 

%% download Single cycle relative times and the folder 'Functional_Data' in
% the current folder. Data are automatically taken from these files.

%% variables.

Current_folder = pwd; % identify current folder.

prompt = {'Number of POIs','Dwell time (ms)','Fill time (ms)', 'Number of trials','green or red'};
dlg_title = 'POIs exp variables';
num_lines = 1;
defaultans = {'','0.004','0.0245','3','green'};
ans_POIs = inputdlg(prompt,dlg_title,num_lines,defaultans);

Numb_POIs = str2num(ans_POIs{1,1});
Dwell_time = str2num(ans_POIs{2,1}); % in ms.
Fill_time = str2num(ans_POIs{3,1}); % in ms.
Numb_trials = str2num(ans_POIs{4,1});
Channel = ans_POIs{5,1};
Smoothing_fact = 10;
POI_coordinates = importdata('ROI.dat');
POI_coordinates.POIs_X = POI_coordinates.data(:,5);
POI_coordinates.POIs_Y = POI_coordinates.data(:,6);
POI_coordinates.POIs_Z = POI_coordinates.data(:,7);

%% import POI data and regroup values per POI per trial and per POI.

fdpath = 'Functional_Data\'; % import POIs data from this folder.
clear POIs_data_trial POIs_data

if strcmp(Channel, 'green') == 1
    channel_let = 'B';
elseif strcmp(Channdl,'red') == 1
    channel_let = 'A';
    
end

% POI_data_trial: 1 x Number of PoI, {i} = Number of cycles x number of trials
% POI_data: (Number of cycles*Number of trials) x Number of PoI 
%           Trial 1;
%             ...
%           Trial N
% Filename format: ROI- $ROI-Number as NNN$ _POI_Dwell-4.0us_Ch 
%                  $A(Red) or B(Green)$ _Trial- $Trial Number as NN$
% For loops not needed - used sprintf instead


for poi_n = 1:Numb_POIs
    for trial_n = 1:Numb_trials
        POIs_data_trial{poi_n}(:,trial_n) = importdata([fdpath,'ROI-',sprintf('%03d',poi_n),... 
            '_POI_Dwell-4.0us_Ch', channel_let,'_Trial-',sprintf('%02d',trial_n),'.dat']);
        POIs_data{trial_n}(:,poi_n) = POIs_data_trial{poi_n}(:,trial_n);
    end
end
POIs_data = vertcat(POIs_data{:}); % every trials are regroup per point.

% POI_data_trial: 1 x Number of PoI, {i} = Number of cycles x number of trials
% POI_data: (Number of cycles*Number of trials) x Number of PoI 
%           Trial 1;
%             ...
%           Trial N

%% time stamp for the experiment.

Numb_cycle = size(POIs_data,1)/Numb_trials; % Alternatively read from header file.
POIs_data_time = zeros(size(POIs_data));
POIs_cycletime_trial1 = zeros(Numb_cycle, Numb_POIs);

POIs_cycletime = importdata('Single cycle relative times.txt'); % values are in microseconds.

if exist('Time_between_trial', 'var')

    if Numb_POIs == size(POIs_cycletime, 1)
        cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
        acquisition_rate = 1000/cycle_length; % in Hz.
        % determine the time stamps for each points for one trial.
        for i = 1:Numb_POIs
            for j = 1:Numb_cycle
                POIs_cycletime_trial1(j,i) = (POIs_cycletime(i,1)/1000) + (cycle_length * (j-1));
            end
        end
        % and for all the trials.
        for i = 1:Numb_POIs
            for j = 1:Numb_trials
                if j == 1
                    POIs_data_time(1:Numb_cycle,i) = POIs_cycletime_trial1(:,i);
                end
                if j > 1
                    POIs_data_time((Numb_cycle*(j-1))+1:Numb_cycle*j,i) = POIs_data_time(Numb_cycle*(j-1),i)+Time_between_trial(1,j-1)+POIs_cycletime_trial1(:,i);
                end
            end
        end
        TimeTotal_imaging = POIs_data_time(end,end);
    end

    if Numb_POIs < size(POIs_cycletime, 1)
        cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
        Numb_repetition = size(POIs_cycletime, 1)/Numb_POIs; % Number of cycles without movement correction. triggered.
        Numb_cycle_MC = round(Numb_cycle/Numb_repetition); % number of cycle with movement correction.
        % time stamps for the first trial.
        for i = 1:Numb_POIs
            for ii = 1:Numb_cycle_MC
                for j = 1:Numb_repetition
                    POIs_cycletime_trial1(j+((Numb_repetition*ii-Numb_repetition)),i) =(POIs_cycletime(i+(Numb_POIs*j-Numb_POIs),1)/1000) + (cycle_length*ii-cycle_length); % in ms.
                end
            end
        end
        if size(POIs_cycletime_trial1, 1) > size(POIs_data, 1)/Numb_trials; % Numb_cycle_MC is an approximate so extra data points could be calculated. The matrix is reshape to fit the real number of data values.
            POIs_cycletime_trial1 = POIs_cycletime_trial1(1:size(POIs_data,1)/Numb_trials,:);
        end
        acquisition_rate = (1000*size(POIs_cycletime_trial1,1))/POIs_cycletime_trial1(end,end); % in Hz.
        % and for all the trials.
        POIs_data_time(1:size(POIs_cycletime_trial1,1),:) = POIs_cycletime_trial1(:,:);
        for i = 1:Numb_POIs;
            for ii = 2:Numb_trials;
                POIs_data_time((size(POIs_cycletime_trial1,1)*(ii-1))+1:size(POIs_cycletime_trial1,1)*ii,i) = POIs_data_time(size(POIs_cycletime_trial1,1)*(ii-1),i)+Time_between_trial(1,ii-1)+POIs_cycletime_trial1(:,i);
            end
        end
        TimeTotal_imaging = POIs_data_time(end,end); % in ms.
    end

end

if exist('Time_between_trial') == 0
    if Numb_POIs == size(POIs_cycletime, 1);
        cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
        acquisition_rate = 1000/cycle_length; % in Hz.
        TimeTotal_imaging = ((cycle_length*Numb_cycle)*Numb_trials); % in ms.
        for i = 1:Numb_POIs;
            for j = 1:size(POIs_data_time,1);
                POIs_data_time(j,i) = (POIs_cycletime(i,1)/1000) + (cycle_length * (j-1));
            end
        end
        TimeTotal_imaging = POIs_data_time(end,end); % in ms.
    end
    if Numb_POIs < size(POIs_cycletime, 1);
        cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
        Numb_repetition = size(POIs_cycletime, 1)/Numb_POIs; % Number of cycles without movement correction. triggered.
        Numb_cycle_MC = round(Numb_cycle/Numb_repetition); % number of cycle with movement correction.
        for i = 1:Numb_POIs;
            for ii = 1:Numb_cycle_MC;
                for j = 1:Numb_repetition;
                    POIs_cycletime_trial1(j+((Numb_repetition*ii-Numb_repetition)),i) =(POIs_cycletime(i+(Numb_POIs*j-Numb_POIs),1)/1000) + (cycle_length*ii-cycle_length); % in ms.
                end
            end
        end
        if size(POIs_cycletime_trial1, 1) > size(POIs_data, 1)/Numb_trials; % Numb_cycle_MC is an approximate so extra data points could be calculated. The matrix is reshape to fit the real number of data values.
            POIs_cycletime_trial1 = POIs_cycletime_trial1(1:size(POIs_data,1)/Numb_trials,:);
        end
        POIs_data_time(1:size(POIs_cycletime_trial1,1),:) = POIs_cycletime_trial1(:,:);
        for i = 1:Numb_POIs;
            for ii = 2:Numb_trials;
                for j = (size(POIs_cycletime_trial1,1)*(ii-1))+1:size(POIs_cycletime_trial1,1)*ii;
                    POIs_data_time(j,i) = POIs_cycletime_trial1(end,i)*(ii-1) + POIs_cycletime_trial1(j-(size(POIs_cycletime_trial1,1)*(ii-1)),i);
                end
            end
        end
        TimeTotal_imaging = POIs_data_time(end,end); % in ms.
    end

    % Time stamp organized in trial.
    POIs_data_time_trial = {};
    for ii = 1:Numb_POIs;
        POIs_data_time_trial{1,ii} = reshape(POIs_data_time(:,ii), [size(POIs_data,1)/Numb_trials,Numb_trials]);
    end
end
cd (Current_folder); % return to main folder for saving analysis.


%% dFF on continuous POIs data.

for i = 1:size(POIs_data, 2); % determine the median grey value for each trace by fitting a kernel function. Can be visualized by using the function histfit and a kernel.
    pd = fitdist(POIs_data(:,i), 'kernel');
    baseline_median(1,i) = median(pd);
    Std_POIs{i} = std(pd);
end
baseline_intensity = repmat(baseline_median,size(POIs_data,1),1); % adjust the number of rows in the matrix to match POIs
dF_POIs_data = POIs_data - baseline_intensity;
dFF_POIs_data = dF_POIs_data./baseline_intensity;
data_smooth = smooth(dFF_POIs_data, Smoothing_fact);
dFF_POIs_data_smooth = reshape(data_smooth, size(dFF_POIs_data,1), size(dFF_POIs_data,2));

% make a figure of all POIs.

y_offset = 1;
figure; hold on
for i = 1:Numb_POIs;
    Y_offset = (y_offset * i) - y_offset;
    plot(POIs_data_time(:,1)./1000, dFF_POIs_data_smooth(:,i) + Y_offset);
end
axis([0 POIs_data_time(end,end)/1000 0 inf]);
xlabel('Time (s)')
title('dFF smooth all POIs');
saveas(gcf, 'dFF smooth all POIs.fig');
hold off


%% dFF for each POIs organized per trial.

for i = 1:Numb_POIs;
    for ii = 1:size(POIs_data_trial{1,i},2);
        pd = fitdist(POIs_data_trial{1,i}(:,ii), 'kernel');
        baseline_median_trial{1,i}(1,ii) = median(pd);
        Std_POIs_trial{1,i}(1,ii) = std(pd);
    end
    baseline_intensity_trial{1,i} = repmat(baseline_median_trial{1,i}, size(POIs_data_trial{1,i},1), 1); % adjust the number of rows in the matrix to match POIs
    dF_POIs_data_trial{1,i} = POIs_data_trial{1,i} - baseline_intensity_trial{1,i};
    dFF_POIs_data_trial{1,i} = dF_POIs_data_trial{1,i}./baseline_intensity_trial{1,i};
    data_smooth = smooth(dFF_POIs_data_trial{1,i}, Smoothing_fact);
    dFF_POIs_data_trial_smooth{1,i} = reshape(data_smooth, size(dFF_POIs_data_trial{1,i}, 1), size(dFF_POIs_data_trial{1,i}, 2));
end
   
%% --------------------------------------------------------------

%% Commented out as channel option buiklt into importdata args
% 
% if strcmp(Channel, 'red') == 1;
%     % regoup each POI per trial.
%     for ii = 1:Numb_POIs;
%         for i = 1:Numb_trials;
%             if ii < 10;
%                 if (i < 10);
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-00',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-00',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%             if ii >= 10 && ii < 100;
%                 if (i < 10);
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-0',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-0',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%             if ii >= 100;
%                 if (i < 10);
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data_trial{ii}(:,i) = importdata(['ROI-',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%         end
%     end
%     
%     % regroup every trial per POI.
%     for i = 1:Numb_trials;
%         for ii = 1:Numb_POIs;
%             if ii < 10;
%                 if (i < 10);
%                     POIs_data{i}(:,ii) = importdata(['ROI-00',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data{i}(:,ii) = importdata(['ROI-00',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%             if ii >= 10 && ii < 100;
%                 if (i < 10);
%                     POIs_data{i}(:,ii) = importdata(['ROI-0',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data{i}(:,ii) = importdata(['ROI-0',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%             if ii >= 100;
%                 if (i < 10);
%                     POIs_data{i}(:,ii) = importdata(['ROI-',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-0',num2str(i),'.dat']);
%                 else
%                     POIs_data{i}(:,ii) = importdata(['ROI-',num2str(ii),'_POI_Dwell-4.0us_ChA_Trial-',num2str(i),'.dat']);
%                 end
%             end
%         end
%     end
%     
%     % time stamp for the experiment.
%     POIs_data = vertcat(POIs_data{:}); % every trials are regroup per point.
%     POIs_data_time = zeros(size(POIs_data,1),size(POIs_data, 2));
%     POIs_cycletime_trial1 = zeros(Numb_cycle, size(POIs_data, 2));
%     
%     if exist('Time_between_trial') == 1;
%         
%         if Numb_POIs == size(POIs_cycletime, 1);
%             cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
%             acquisition_rate = 1000/cycle_length; % in Hz.
%             % determine the time stamps for each points for one trial.
%             for i = 1:Numb_POIs;
%                 for j = 1:Numb_cycle;
%                     POIs_cycletime_trial1(j,i) = (POIs_cycletime(i,1)/1000) + (cycle_length * (j-1));
%                 end
%             end
%             % and for all the trials.
%             for i = 1:Numb_POIs;
%                 for j = 1:Numb_trials;
%                     if j == 1;
%                         POIs_data_time(1:Numb_cycle,i) = POIs_cycletime_trial1(:,i);
%                     end
%                     if j > 1;
%                         POIs_data_time((Numb_cycle*(j-1))+1:Numb_cycle*j,i) = POIs_data_time(Numb_cycle*(j-1),i)+Time_between_trial(1,j-1)+POIs_cycletime_trial1(:,i);
%                     end
%                 end
%             end
%             TimeTotal_imaging = POIs_data_time(end,end);
%         end
%         
%         if Numb_POIs < size(POIs_cycletime, 1);
%             cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
%             Numb_repetition = size(POIs_cycletime, 1)/Numb_POIs; % Number of cycles without movement correction. triggered.
%             Numb_cycle_MC = round(Numb_cycle/Numb_repetition); % number of cycle with movement correction.
%             % time stamps for the first trial.
%             for i = 1:Numb_POIs;
%                 for ii = 1:Numb_cycle_MC;
%                     for j = 1:Numb_repetition;
%                         POIs_cycletime_trial1(j+((Numb_repetition*ii-Numb_repetition)),i) =(POIs_cycletime(i+(Numb_POIs*j-Numb_POIs),1)/1000) + (cycle_length*ii-cycle_length); % in ms.
%                     end
%                 end
%             end
%             if size(POIs_cycletime_trial1, 1) > size(POIs_data, 1)/Numb_trials; % Numb_cycle_MC is an approximate so extra data points could be calculated. The matrix is reshape to fit the real number of data values.
%                 POIs_cycletime_trial1 = POIs_cycletime_trial1(1:size(POIs_data,1)/Numb_trials,:);
%             end
%             acquisition_rate = (1000*size(POIs_cycletime_trial1,1))/POIs_cycletime_trial1(end,end); % in Hz.
%             % and for all the trials.
%             POIs_data_time(1:size(POIs_cycletime_trial1,1),:) = POIs_cycletime_trial1(:,:);
%             for i = 1:Numb_POIs;
%                 for ii = 2:Numb_trials;
%                     POIs_data_time((size(POIs_cycletime_trial1,1)*(ii-1))+1:size(POIs_cycletime_trial1,1)*ii,i) = POIs_data_time(size(POIs_cycletime_trial1,1)*(ii-1),i)+Time_between_trial(1,ii-1)+POIs_cycletime_trial1(:,i);
%                 end
%             end
%             TimeTotal_imaging = POIs_data_time(end,end); % in ms.
%         end
%         
%     end
%     
%     if exist('Time_between_trial') == 0
%         if Numb_POIs == size(POIs_cycletime, 1);
%             cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
%             acquisition_rate = 1000/cycle_length; % in Hz.
%             TimeTotal_imaging = ((cycle_length*Numb_cycle)*Numb_trials); % in ms.
%             for i = 1:Numb_POIs;
%                 for j = 1:size(POIs_data_time,1);
%                     POIs_data_time(j,i) = (POIs_cycletime(i,1)/1000) + (cycle_length * (j-1));
%                 end
%             end
%             TimeTotal_imaging = POIs_data_time(end,end); % in ms.
%         end
%         if Numb_POIs < size(POIs_cycletime, 1);
%             cycle_length = (POIs_cycletime(1,2))/1000; % in ms.
%             Numb_repetition = size(POIs_cycletime, 1)/Numb_POIs; % Number of cycles without movement correction. triggered.
%             Numb_cycle_MC = round(Numb_cycle/Numb_repetition); % number of cycle with movement correction.
%             for i = 1:Numb_POIs;
%                 for ii = 1:Numb_cycle_MC;
%                     for j = 1:Numb_repetition;
%                         POIs_cycletime_trial1(j+((Numb_repetition*ii-Numb_repetition)),i) =(POIs_cycletime(i+(Numb_POIs*j-Numb_POIs),1)/1000) + (cycle_length*ii-cycle_length); % in ms.
%                     end
%                 end
%             end
%             if size(POIs_cycletime_trial1, 1) > size(POIs_data, 1)/Numb_trials; % Numb_cycle_MC is an approximate so extra data points could be calculated. The matrix is reshape to fit the real number of data values.
%                 POIs_cycletime_trial1 = POIs_cycletime_trial1(1:size(POIs_data,1)/Numb_trials,:);
%             end
%             POIs_data_time(1:size(POIs_cycletime_trial1,1),:) = POIs_cycletime_trial1(:,:);
%             for i = 1:Numb_POIs;
%                 for ii = 2:Numb_trials;
%                     for j = (size(POIs_cycletime_trial1,1)*(ii-1))+1:size(POIs_cycletime_trial1,1)*ii;
%                         POIs_data_time(j,i) = POIs_cycletime_trial1(end,i)*(ii-1) + POIs_cycletime_trial1(j-(size(POIs_cycletime_trial1,1)*(ii-1)),i);
%                     end
%                 end
%             end
%             TimeTotal_imaging = POIs_data_time(end,end); % in ms.
%         end
%         % Time stamp organized in trial.
%         POIs_data_time_trial = {};
%         for ii = 1:Numb_POIs;
%             POIs_data_time_trial{1,ii} = reshape(POIs_data_time(:,ii), [size(POIs_data,1)/Numb_trials,Numb_trials]);
%         end
%     end
%     cd (Current_folder); % return to main folder for saving analysis.
%     
%     
%     %% dFF on continuous POIs data.
%     
%     for i = 1:size(POIs_data, 2); % determine the median grey value for each trace by fitting a kernel function. Can be visualized by using the function histfit and a kernel.
%         pd = fitdist(POIs_data(:,i), 'kernel');
%         baseline_median(1,i) = median(pd);
%         Std_POIs{i} = std(pd);
%     end
%     baseline_intensity = repmat(baseline_median,size(POIs_data,1),1); % adjust the number of rows in the matrix to match POIs
%     dF_POIs_data = POIs_data - baseline_intensity;
%     dFF_POIs_data = dF_POIs_data./baseline_intensity;
%     data_smooth = smooth(dFF_POIs_data, Smoothing_fact);
%     dFF_POIs_data_smooth = reshape(data_smooth, size(dFF_POIs_data,1), size(dFF_POIs_data,2));
%     
%     % make a figure of all POIs.
%     
%     y_offset = 1;
%     figure; hold on
%     for i = 1:Numb_POIs;
%         Y_offset = (y_offset * i) - y_offset;
%         plot(POIs_data_time(:,1)./1000, dFF_POIs_data_smooth(:,i) + Y_offset);
%     end
%     axis([0 POIs_data_time(end,end)/1000 0 inf]);
%     xlabel('Time (s)')
%     title('dFF smooth all POIs');
%     saveas(gcf, 'dFF smooth all POIs.fig');
%     hold off
%     
%     
%     %% dFF for each POIs organized per trial.
%     
%     for i = 1:Numb_POIs;
%         for ii = 1:size(POIs_data_trial{1,i},2);
%             pd = fitdist(POIs_data_trial{1,i}(:,ii), 'kernel');
%             baseline_median_trial{1,i}(1,ii) = median(pd);
%             Std_POIs_trial{1,i}(1,ii) = std(pd);
%         end
%         baseline_intensity_trial{1,i} = repmat(baseline_median_trial{1,i}, size(POIs_data_trial{1,i},1), 1); % adjust the number of rows in the matrix to match POIs
%         dF_POIs_data_trial{1,i} = POIs_data_trial{1,i} - baseline_intensity_trial{1,i};
%         dFF_POIs_data_trial{1,i} = dF_POIs_data_trial{1,i}./baseline_intensity_trial{1,i};
%         data_smooth = smooth(dFF_POIs_data_trial{1,i}, Smoothing_fact);
%         dFF_POIs_data_trial_smooth{1,i} = reshape(data_smooth, size(dFF_POIs_data_trial{1,i}, 1), size(dFF_POIs_data_trial{1,i}, 2));
%     end
%end
