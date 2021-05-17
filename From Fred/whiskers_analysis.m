%% import measurements and whiskers data for the whiskers tracking analysis (openwiki, Janelia farm).

currentfolder = cd;
prompt = {'measurements file name','whiskers file name'};
dlg_title = 'Whiskers tracking file name';
num_lines = 1;
defaultans = {'',''};
ans_filename = inputdlg(prompt,dlg_title,num_lines,defaultans);
Measure_path = fullfile(currentfolder, ans_filename{1,1});
Whiskers_path = fullfile(currentfolder, ans_filename{2,1});

Whiskers.measurements = LoadMeasurements(Measure_path);
%[Whiskers.whiskers, Whiskers.format] = LoadWhiskers(Whiskers_path);

clear ans_filename currentfolder defaultans dlg_title Measure_path Whiskers_path num_lines prompt

%% sort whiskers angle in function of their label (-1 are not considered as whiskers).

Whiskers.label = [Whiskers.measurements.label].';
Whiskers.frameid = [Whiskers.measurements.fid].';
Whiskers.angle = [Whiskers.measurements.angle].';

% number of whiskers id.
Whiskers.whiskers_numb = unique(Whiskers.label);

for id = 0:max(Whiskers.whiskers_numb);
    [row, ~] = find(Whiskers.label == id);
    Whiskers.row_label{id+1} = row; % whisker id start from 1.
    Whiskers_frame = Whiskers.frameid(row,1);
    Whiskers_angle = Whiskers.angle(row,1);
    Whiskers_angle_smooth = smooth(Whiskers_angle, 10);
    Whiskers_time = zeros(size(Whiskers_frame,1),1);
    for i = 1:size(Whiskers_frame,1);
        [row, ~] = find(Whiskers.timestamp(:,1) == Whiskers_frame(i,1)+1);
        Whiskers_time(i,1) = Whiskers.timestamp(row,2);
    end 
    assignin('base', sprintf('Whiskers_frame_%d', id), Whiskers_frame);
    assignin('base', sprintf('Whiskers_angle_%d', id), Whiskers_angle);
    assignin('base', sprintf('Whiskers_angle_smooth_%d', id), Whiskers_angle_smooth);
    assignin('base', sprintf('Whiskers_time_%d', id), Whiskers_time);
    clear Whiskers_frame Whiskers_angle Whiskers_time Whiskers_angle_smooth row i id
end

Whiskers = rmfield(Whiskers, 'measurements'); % remove the whiki file measurements to save space.







