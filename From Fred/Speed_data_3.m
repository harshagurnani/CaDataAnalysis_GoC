% updated 02/12/2016

%% import speed encoder data from rig3.
% in the speed data file, time is in microseconds and speed in rpm (1 rotation =
% 50 cm). Import speed data from different modes of acquisition and
% determine the time trials from the speed data and the imaging acquisition
% time for one trial. This script needs to be run first in order to
% determine the time between trials and the correct time stamps for the
% imaging data. 

prompt = {'TS, POIs or patch'};
dlg_title = 'Acquisition mode';
num_lines = 1;
defaultans = {''};
ans_mode = inputdlg(prompt,dlg_title,num_lines,defaultans);
Mode = ans_mode{1,1};

if strcmp(Mode, 'TS') == 1;
    Speed_data = importdata('Speed data.txt');
    Speed = Speed_data.data(:,[1,2]);
else
    Speed_data = importdata('Speed data 001.txt');
    Speed = Speed_data.data(:,[1,2]);
end
% exp parameters for patch experiments.
if strcmp(Mode, 'patch') == 1;
    Cycle = importdata('Single cycle relative times.txt');
    Time_cycle = Cycle(1,2); % in microseconds.
    prompt = {'Number of cycles'};
    dlg_title = 'Patch parameters';
    num_lines = 1;
    defaultans = {'',};
    ans_patch = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Numb_cycle = str2num(ans_patch{1,1});
    Time_imaging_per_trial = (Time_cycle*Numb_cycle)/1000; % in ms.
end
% exp parameters for POIs experiments.
if strcmp(Mode, 'POIs') == 1;
    prompt = {'Number of cycles', 'Number of POIs'};
    dlg_title = 'POIs parameters';
    num_lines = 1;
    defaultans = {'',''};
    ans_POIs = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Cycle = importdata('Single cycle relative times.txt');
    if size(Cycle(:,1),1) == str2num(ans_POIs{2,1});
        Time_cycle = Cycle(1,2); % in microseconds.
        Numb_cycle = str2num(ans_POIs{1,1});
        Time_imaging_per_trial = (Time_cycle*Numb_cycle)/1000; % in ms.
    else
        cycle_length = Cycle(1,2);
        Numb_repetition = size(Cycle(:,1),1)/str2num(ans_POIs{2,1});
        Numb_cycle_MC = round(str2num(ans_POIs{1,1})/Numb_repetition);
        Time_imaging_per_trial = (cycle_length*Numb_cycle_MC)/1000; % in ms.
    end
end

prompt = {'number of trials','speed threshold (cm.s)'};
dlg_title = 'number of trials ?';
num_lines = 1;
defaultans = {'','2.5'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

Numb_trials = str2num(answer{1,1});
Speed_threshold = str2num(answer{2,1}); % in cm/s. threshold for locomotion activity.
clear SpeedTimeMatrix row counter FlagRow SpeedDataMatrix
Speed(:,2) = Speed(:,2)*-1;

%% Speed data for one trial.
row = 0;
if Numb_trials == 1;
    for row = row+1:length(Speed)-1; % find the end of the trial.
        if Speed(row, 1) > Speed(row+1, 1);
            SpeedTime_temp = Speed(1:row, 1);
            SpeedData_temp = Speed(1:row, 2);
            SpeedTimeMatrix = SpeedTime_temp; % create a matrix with the time in microseconds.
            SpeedDataMatrix = SpeedData_temp;
        end
    end
    if Speed(row, 1) < Speed(row+1, 1);
        SpeedTimeMatrix = Speed(:,1);
        SpeedDataMatrix = Speed(:,2);
    end
end


%% Identify the trials and the time between trials.
if Numb_trials > 1;
    %clear SpeedTimeTrial SpeedDataTrial SpeedTimeTrial_temp SpeedDataTrial_temp
    counter = 0;
    row = 0;
    for row = row+1:length(Speed)-1; % find the different trials in speed data.
        if Speed(row, 1) > Speed(row+1, 1);
            counter = counter + 1;
            if counter == 1;
                SpeedTimeTrial_temp_all{counter} = Speed(1:row, 1);
                a = row;
            end
            if counter > 1;
                SpeedTimeTrial_temp_all{counter} = Speed(a+1:row,1);
                a = row;
            end
        end
    end
    SpeedTimeTrial_temp_all{counter+1} = Speed(a+1:end,1);
    if size(SpeedTimeTrial_temp_all,2) < Numb_trials*2;
        SpeedTimeTrial_temp_all{counter+2} = 0;
    end
    clear a
    
    %% add up the trials and the time between trials for rig 3 speed data update.
    SpeedTimeTrial_temp = {};
    
    if size(SpeedTimeTrial_temp,2) == Numb_trials; % correspond to the older version of speed data on rig3.
        SpeedTimeTrial_temp = SpeedTimeTrial_temp_all;
    else
        size(SpeedTimeTrial_temp_all,2) > Numb_trials; % add the pause between trial.
        ii = 2:2:size(SpeedTimeTrial_temp_all,2);
        i = 0;
        for j = 1:2:size(SpeedTimeTrial_temp_all,2);
            i = i+1;
            SpeedTimeTrial_temp{1,i} = zeros(size(SpeedTimeTrial_temp_all{1,j},1)+size(SpeedTimeTrial_temp_all{1,ii(i)},1),1);
            SpeedTimeTrial_temp{1,i}(1:size(SpeedTimeTrial_temp_all{1,j},1),1) = SpeedTimeTrial_temp_all{1,j};
            SpeedTimeTrial_temp{1,i}(size(SpeedTimeTrial_temp_all{1,j},1)+1:end) = SpeedTimeTrial_temp_all{1,j}(end)+SpeedTimeTrial_temp_all{1,ii(i)};
        end
    end
    clear i ii j
    
    
    
    %% calculate the time between trials in function of the imaging time for one trial.
    
    Time_between_trial = zeros(1,size(SpeedTimeTrial_temp,2));
    for i = 1:size(SpeedTimeTrial_temp,2);
        Time_between_trial(1,i) = SpeedTimeTrial_temp{1,i}(end)/1000 - Time_imaging_per_trial; % in ms
    end
    
    %% Create matrices for the speed data and the speed time stamps.
    SpeedTimeTrial = SpeedTimeTrial_temp;
    SpeedTimeMatrix = zeros(length(Speed),1); % create a matrix with the time in microseconds.
    SpeedDataMatrix = Speed(:,2);
    for i = 1:Numb_trials;
        if i == 1;
            SpeedTimeMatrix(1:length(SpeedTimeTrial{1,i})) = SpeedTimeTrial{1,i}(:,1);
            counter = length(SpeedTimeTrial{1,i});
        end
        if i > 1;
            SpeedTimeMatrix(counter+1:counter+length(SpeedTimeTrial{1,i}),1) = SpeedTimeMatrix(counter,1)+SpeedTimeTrial{1,i}(:,1);
            counter = counter + length(SpeedTimeTrial{1,i});
        end
    end
end

if size(SpeedDataMatrix,1) < size(SpeedTimeMatrix,1);
    SpeedTimeMatrix = SpeedTimeMatrix(1:size(SpeedDataMatrix,1),1);
end

TimeTotal_speed = SpeedTimeMatrix(end)/1000; % in ms.


%% Define period of rest, activity and locomotion.
SpeedTimeMatrix = SpeedTimeMatrix/1000; % in ms.
SpeedDataMatrix = (SpeedDataMatrix(:,1)*50)/60; % in cm/s.
SpeedDataMatrix_smooth = smooth(SpeedDataMatrix, 10);

% time series of speed data

Locomotion_data = timeseries(SpeedDataMatrix_smooth, SpeedTimeMatrix); % change x axis.
set(Locomotion_data, 'name', 'Locomotion data');
plot(Locomotion_data);
xlabel('Time (milliseconds)');
ylabel('cm/s');
axis([0 inf -5 inf]);
savefig('Locomotion data.fig');

% define locomotion period.

clear rloc Time_loc Time_loc_range Numb_loc_periods
[rloc, ~] = find(SpeedDataMatrix > Speed_threshold); % find points above a threshold corresponding to locomotion in cm/s.
a = isempty(rloc);
if a == 1;
    Numb_loc_periods = 0;
end
if a == 0;
    Time_loc = zeros(size(rloc)); % create a matrix with the time stamp of these values in ms.
    for i = 1:size(Time_loc);
        Time_loc(i, 1) = SpeedTimeMatrix(rloc(i, 1),1);
    end
    counter = 1;
    for j = 1:size(Time_loc, 1)-1; % determine the range in ms of locomotion periods.
        if Time_loc(j+1) > Time_loc(j) + 2000; % locomotion periods spaced by 2000 ms.
            if counter == 1;
                Time_loc_range{counter}(1,1) = Time_loc(1);
                Time_loc_range{counter}(1,2) = Time_loc(j);
                jj = Time_loc(j+1);
            else
                Time_loc_range{counter}(1,1) = jj;
                Time_loc_range{counter}(1,2) = Time_loc(j);
                jj = Time_loc(j+1);
            end
            counter = counter+1;
        end
    end
    if counter > 1;
        Time_loc_range{counter}(1,1) = jj; % find the last locomotion period
        Time_loc_range{counter}(1,2) = Time_loc(end,1);
    end
    b = exist('Time_loc_range'); % if only one locomotion period.
    if b == 0;
        Time_loc_range{counter}(1,1) = Time_loc(1, 1);
        Time_loc_range{counter}(1,2) = Time_loc(end, 1);
    end
    Numb_loc_periods = size(Time_loc_range, 2);
    for k = 1:size(Time_loc_range, 2);
        Loc_dur{k} = (Time_loc_range{1,k}(1,2) - Time_loc_range{1,k}(1,1)) / 1000; % in s.
    end
    Loc_duration = zeros(size(Time_loc_range));
    for k = 1:size(Time_loc_range, 2);
        Loc_duration(1,k) = Loc_dur{1,k}(1,1);
    end
    Mean_locomotion_period = mean(Loc_duration,2); % in s.
end
clear a b


% dynamic bar graph for the speed data.

prompt = {'yes=1, no=0'};
dlg_title = 'Dynamic bar graph?';
num_lines = 1;
defaultans = {'0'};
ans_dynbargraph = inputdlg(prompt,dlg_title,num_lines,defaultans);
Flag_dynbargraph = str2num(ans_dynbargraph{1,1});

if Flag_dynbargraph == 1;
    prompt = {'Step size', 'Speed up factor'};
    dlg_title = 'Video parameters';
    num_lines = 1;
    defaultans = {'1', '1'};
    ans_stepsize = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Speed_acq_rate = (size(SpeedDataMatrix,1)*1000)/TimeTotal_speed; % in Hz.
    videorate = (Speed_acq_rate/str2num(ans_stepsize{1,1}))*str2num(ans_stepsize{2,1});
    %[dynamicbargraph] = Dynamicbargraph(Data, framerate, stepsize, Ylabel)
    [dynamicbargraph] = Dynamicbargraph(SpeedDataMatrix_smooth, round(videorate), str2num(ans_stepsize{1,1}), 'cm/s');
end



