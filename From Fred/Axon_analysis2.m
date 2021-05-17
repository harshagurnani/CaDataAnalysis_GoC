%% create figures with dFF, traces, in parallel with locomotion and puff
% experiments.
% Determine bursts time after a puff.

prompt = {'event detection?','Puff?'};
dlg_title = 'Flags; yes=1, no=0';
num_lines = 1;
defaultans = {'0',''};
ans_flags = inputdlg(prompt,dlg_title,num_lines,defaultans);

Flag_eventdetection = str2num(ans_flags{1,1}); % if flag = 1, event detection is performed.
Flag_puff = str2num(ans_flags{2,1}); % if = 0, no air puff experiment.
Loco_data = exist('Locomotion_data');

%% create timeseries for each axon data, figure and smooth data.

Numb_total_axon = size(Axon_dFF,2);
for i = 1:Numb_total_axon; % for each axon.
    Axon_ts{i} = timeseries(Axon_dFF(:,i), TimeAxon(:,i)); % time in ms.
end

s = 10; % smoothing factor.
[m, n] = size(Axon_dFF);
Data_smooth = smooth(Axon_dFF, s);
Axon_dFF_smooth = reshape(Data_smooth,m,n);

% figure with traces for each axon.

y_offset = 0.5;
figure; hold on
for i = 1:size(Axon_dFF_smooth, 2);
    Y_offset = (y_offset * i) - y_offset;
    plot(TimeAxon(:,1)./1000, Axon_dFF_smooth(:,i) + Y_offset);
end
axis([0 inf -2 inf]);
xlabel('Time (s)');
title('dFF smooth all axons');
saveas(gcf, 'dFF smooth all axons_2.fig');
hold off

% figure with the dFF overtime and the locomotion period.
figure; hold on
subplot(2,1,1);
Axon_dFF_smooth_norm = bsxfun(@rdivide, Axon_dFF_smooth, max(Axon_dFF_smooth));
imagesc(Axon_dFF_smooth_norm', [0 1]);
ax = gca;
ax.XAxis.Visible = 'off';
colormap(jet);
if size(Axon_dFF_smooth_norm, 2)>1;
    axis([0 size(Axon_dFF_smooth_norm, 1) 1 size(Axon_dFF_smooth_norm, 2)]);
end
axis ij;
if Loco_data == 1;
    subplot(2,1,2);
    plot(SpeedTimeMatrix(:,1)./1000, SpeedDataMatrix_smooth(:,1), 'k', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('cm/s');
    axis([0 inf -5 inf])
end
saveas(gcf, 'PFs Locomotion.fig');
hold off
clear ax

%% change of dFF for each axon during each locomotion period.
if Loco_data == 1;
    if Numb_loc_periods > 0;
        
        y_offset = 0.5;
        figure; hold on
        for i = 1:size(Axon_dFF_smooth, 2);
            Y_offset = (y_offset * i) - y_offset;
            plot(TimeAxon(:,1)./1000, Axon_dFF_smooth(:,i) + Y_offset);
        end
        axis([0 inf -2 inf]);
        for i = 1:Numb_loc_periods;
            p = patch([Time_loc_range{1,i}(1,1)/1000 Time_loc_range{1,i}(1,2)/1000 Time_loc_range{1,i}(1,2)/1000 Time_loc_range{1,i}(1,1)/1000], [-2 -2 y_offset*Numb_total_axon+1 y_offset*Numb_total_axon+1], 'k');
            set(p,'FaceAlpha',0.2);
            set(p, 'EdgeColor', 'none');
        end
        %ax = gca;
        %ax.YAxis.Visible = 'off';
        xlabel('Time (s)');
        title('dFF smooth all axons');
        saveas(gcf, 'dFF smooth all axons and loc periods.fig');
        hold off
        
        clear row Row dFF_smooth_loc Axon_dFF_loc Axon_dFF_smooth_loc Loc_dim Axon_dFF_loc Axon_dFF_smooth_loc dFF_smooth_loc
        
        baseline_loc = 3000; % in ms, time before and after a locomotion period.
        
        for i = 1:Numb_loc_periods;
            for j = 1:Numb_total_axon;
                [row, ~] = find(TimeAxon(:,j) > (Time_loc_range{1,i}(1,1))-baseline_loc & TimeAxon(:,j) < (Time_loc_range{1,i}(1,2))+baseline_loc);
                Row{i,j} = row;
                Axon_dFF_loc{i,j} = Axon_dFF(Row{i,j},j);
                Axon_dFF_smooth_loc{i,j} = Axon_dFF_smooth(Row{i,j},j);
            end
        end
        Loc_dim = zeros(size(Axon_dFF_loc)); % all the traces are reshape to have the same length.
        for i = 1:Numb_loc_periods;
            for j = 1:Numb_total_axon;
                Loc_dim(i,j) = size(Axon_dFF_loc{i,j},1);
            end
        end
        for i = 1:Numb_loc_periods;
            M = min(Loc_dim(i,:));
            for j = 1:Numb_total_axon;
                Axon_dFF_loc{i,j} = Axon_dFF_loc{i,j}(1:M,1);
                Axon_dFF_smooth_loc{i,j} = Axon_dFF_smooth_loc{i,j}(1:M,1);
            end
        end
        
        % make one figure for each locomotion period. norm for the higher dFF
        % per trace.
        for i = 1:Numb_loc_periods;
            figure; hold on
            dFF_smooth_loc = zeros(size(Axon_dFF_smooth_loc{i,1},1),Numb_total_axon);
            for j = 1:Numb_total_axon;
                dFF_smooth_loc(:,j) = Axon_dFF_smooth_loc{i,j}(:,1);
            end
            dFF_smooth_loc_norm = dFF_smooth_loc;
            dFF_smooth_loc_norm = bsxfun(@rdivide, dFF_smooth_loc, max(Axon_dFF_smooth)); % normalized for the max value over the recording.
            subplot(2,1,1);
            imagesc(dFF_smooth_loc_norm', [0 1]);
            colormap(jet);
            ax = gca;
            ax.XAxis.Visible = 'off';
            %axis([0 size(dFF_smooth_loc_norm, 1) 1 size(dFF_smooth_loc_norm, 2)+1])
            axis ij;
            subplot(2,1,2);
            [row, ~] = find(SpeedTimeMatrix(:,1) > (Time_loc_range{1,i}(1,1))-baseline_loc & SpeedTimeMatrix(:,1) < (Time_loc_range{1,i}(1,2))+baseline_loc);
            % 2000 ms before and after the locomotion period.
            plot(SpeedTimeMatrix(row,1)./1000, SpeedDataMatrix_smooth(row, 1), 'k', 'LineWidth', 1);
            axis([SpeedTimeMatrix(row(1,1),1)./1000 inf -5 inf]);
            xlabel('Time (s)');
            saveas(gcf, sprintf('Locomotion%d.fig',i));
            hold off
        end
        clear ax
        
        % plot dFF traces during each locomotion period for all axons.
        
        for i = 1:Numb_loc_periods;
            figure; hold on
            y_offset = 0.1;
            subplot(2,1,1);
            for ii = 1:Numb_total_axon;
                hold on
                Y_offset = (y_offset * ii) - y_offset;
                x = 1:size(Axon_dFF_smooth_loc{i,ii}(:,1),1);
                plot(x', Axon_dFF_smooth_loc{i,ii}(:,1) + Y_offset);
                axis([0 inf -0.2 inf]);
                hold off
            end
            subplot(2,1,2);
            hold on
            [row, ~] = find(SpeedTimeMatrix(:,1) > (Time_loc_range{1,i}(1,1))-baseline_loc & SpeedTimeMatrix(:,1) < (Time_loc_range{1,i}(1,2))+baseline_loc);
            x = (SpeedTimeMatrix(row, 1)./1000)-(SpeedTimeMatrix(row(1,1), 1)./1000);
            plot(x', SpeedDataMatrix_smooth(row, 1), 'k', 'LineWidth', 1);
            axis([0 inf -5 inf]);
            xlabel('Time (s)');
            ylabel('Speed (cm/s)');
            hold off
            saveas(gcf, sprintf('PFs locomotion period_%d.fig',i));
            hold off
            clear x row
        end
    end
    
     % quantify timing of bursts at the beginning of a locomotion period.
     a = exist('Frame_rate');
     if a == 1;
         acquisition_rate = Frame_rate;
     end
     clear a
        % detected if movement below the Speed_threshold in a 2000 ms
        % window before the locomotion period detection. 
        Burstloc.Time_loc_start = cell(1,Numb_loc_periods);
        Burstloc.row_start = cell(1,Numb_loc_periods);
        for i = 1:Numb_loc_periods;
            [row, ~] = find(SpeedTimeMatrix(:,1) <= Time_loc_range{1,i}(1,1) & SpeedTimeMatrix(:,1) > Time_loc_range{1,i}(1,1)-2000); % find the row corresponding to 2000 ms before the locomotion period.
            Burstloc.Speed_abs = abs(SpeedDataMatrix_smooth(row, 1)); % absolute value to catch even if the mouse moved backward.
            [row_low, ~] = find(Burstloc.Speed_abs > 1); % threshold set at 1 cm/s.
            Burstloc.Time_loc_start{1,i} = SpeedTimeMatrix(row(1,1)+row_low(1,1)+1, 1); % time in ms when the wheel moves before a detected locomotion period.
            Burstloc.row_start{1,i} = row(1,1)+row_low(1,1)+1; % row in speed data where the burst analysis will be performed. 
        end
        for i = 1:Numb_loc_periods;
            if Burstloc.Time_loc_start{1,i} < 1000; % if the locomotion period starts less than 1000 ms after the imaging, this period is not analyzed. 
            Burstloc.Time_loc_start{1,i} = NaN;
            end
        end 
     %[burst_time_threshold, Row_event, Row_max, Fit_function, Onset, X, K,...
            %burst_onset_f, burst_onset_w] = burst_detection_trigger(data,...
            %data_time, acquisition_rate, trigger_time, baseline, detection_window);
            %values are in ms.
            Burstloc.Burst_time_threshold = zeros(Numb_trials,Numb_total_axon); % create a matrix that will contain the time of the first burst when it reaches the threshold.
            Burstloc.burst_time_fit = NaN(Numb_loc_periods,Numb_total_axon); % create a matrix that will contain the time of the first burst during the puff for each axons.
            Burstloc.burst_time_w = zeros(Numb_loc_periods,Numb_total_axon);
            Burstloc.row_max = NaN(Numb_loc_periods,Numb_total_axon);
            Burstloc.row_event = zeros(Numb_loc_periods,Numb_total_axon);
            Burstloc.fit_function = cell(Numb_loc_periods,Numb_total_axon);
            Burstloc.onset = zeros(Numb_loc_periods,Numb_total_axon);
            Burstloc.X_f = cell(Numb_loc_periods,Numb_total_axon);
            Burstloc.K_w = cell(Numb_loc_periods,Numb_total_axon);
            
            for i = 1:Numb_loc_periods;
                for ii = 1:Numb_total_axon;
                    if isnan(Burstloc.Time_loc_start{1,i}) == 1;
                        Burstloc.Burst_time_threshold(i,ii) = NaN;
                        Burstloc.burst_time_fit(i,ii) = NaN;
                        Burstloc.burst_time_w(i,ii) = NaN;
                        Burstloc.row_max(i,ii) = NaN;
                        Burstloc.row_event(i,ii) = NaN;
                        Burstloc.fit_function{i,ii} = [];
                        Burstloc.onset(i,ii) = NaN;
                        Burstloc.X_f{i,ii} = NaN;
                        Burstloc.K_w{i,ii} = NaN;
                    else
                        [burst_time_threshold, Row_event, Row_max, Fit_function, Onset, X, K, burst_onset_f, burst_onset_w] = burst_detection_trigger(Axon_dFF_smooth(:,ii), TimeAxon(:,ii), acquisition_rate, Burstloc.Time_loc_start{1,i}, 500, 3000); % values are in ms.
                        Burstloc.Burst_time_threshold(i,ii) = burst_time_threshold;
                        Burstloc.burst_time_fit(i,ii) = burst_onset_f;
                        Burstloc.burst_time_w(i,ii) = burst_onset_w;
                        Burstloc.row_max(i,ii) = Row_max;
                        Burstloc.row_event(i,ii) = Row_event;
                        Burstloc.fit_function{i,ii} = Fit_function;
                        Burstloc.onset(i,ii) = Onset;
                        Burstloc.X_f{i,ii} = X;
                        Burstloc.K_w{i,ii} = K;
                    end
                end
            end
            
            % figures showing where the onset has been taken.
            
    prompt = {'check traces for onset?'};
    dlg_title = 'yes=1, no=0';
    num_lines = 1;
    defaultans = {'1'};
    ans_onset = inputdlg(prompt,dlg_title,num_lines,defaultans);
    ans_onset = str2num(ans_onset{1,1});
    if ans_onset == 1;
        close all
        for i = 1:Numb_total_axon;
            for ii = 1:Numb_loc_periods;
                flag_w = 0;
                flag_fit = 1;
                a = isfield(Burstloc, 'fit_function');
                if a == 1;
                    f = Burstloc.fit_function{ii,i};
                    flag_fit = isempty(f);
                    x = Burstloc.X_f{ii,i};
                end
                b = isfield(Burstloc, 'K_w');
                if b == 1
                    k = Burstloc.K_w{ii,i};
                    if k > 0;
                        flag_w = 1;
                    end
                end
                subplot(2,1,1);
                if flag_fit == 0;
                    plot(f, x', Axon_dFF_smooth(Burstloc.row_event(ii,i):Burstloc.row_max(ii,i),i), 'k');
                    if round(Burstloc.onset(ii, i)) > 0 & round(Burstloc.onset(ii, i)) < size(Axon_dFF_smooth,1);
                        hold on
                        plot(round(Burstloc.onset(ii,i)), Axon_dFF_smooth(Burstloc.row_event(ii,i)+round(Burstloc.onset(ii,i)),i), 'r*');
                        hold off
                    end
                end
                title(['fit detection', sprintf('Axon %d', i), sprintf('Locomotion period %d',ii)]);
                subplot(2,1,2);
                if Burstloc.row_max(ii,i) > 0;
                    plot(Axon_dFF_smooth(Burstloc.row_event(ii,i):Burstloc.row_max(ii,i),i), 'k');
                    if flag_w == 1;
                        hold on
                        plot(k-Burstloc.row_event(ii,i)+1, Axon_dFF_smooth(k,i), 'r*');
                        hold off
                    end
                end
                title('sliding window detection');
                clear f x k flag_fit a b
                waitfor(gcf);
            end
        end
    end
    clear flag_w
    
    % burst time is substracted from the start of the locomotion period.
    for i = 1:Numb_total_axon;
        for ii = 1:Numb_loc_periods;
            Burstloc.burst_time_fit(ii, i) = Burstloc.burst_time_fit(ii, i) - Burstloc.Time_loc_start{1,ii}; % in ms.
        end
    end
    for i = 1:(size(Burstloc.burst_time_fit,1)); % remove operation performed on NaN values.
        for ii = 1:(size(Burstloc.burst_time_fit,2));
            if Burstloc.burst_time_fit(i,ii) == -Burstloc.Time_loc_start{1,i};
                Burstloc.burst_time_fit(i,ii) = NaN;
            end
        end
    end
    
    for i = 1:Numb_total_axon;
        for ii = 1:Numb_loc_periods;
            Burstloc.burst_time_w(ii, i) = Burstloc.burst_time_w(ii, i) - Burstloc.Time_loc_start{1,ii}; % in ms.
        end
    end
    for i = 1:(size(Burstloc.burst_time_w,1)); % remove operation performed on NaN values.
        for ii = 1:(size(Burstloc.burst_time_w,2));
            if Burstloc.burst_time_w(i,ii) == -Burstloc.Time_loc_start{1,i};
                Burstloc.burst_time_w(i,ii) = NaN;
            end
        end
    end
    
    % creates table containing the burst data.
    clear Burstloc.col_var_name Burstloc.row_var_name
    for i = 1:Numb_total_axon;
        Burstloc.col_var_name(1,i) = {sprintf('axon_%d',i)};
    end
    for i = 1:Numb_loc_periods;
        Burstloc.row_var_name(1,i) = {sprintf('locomotion period_%d',i)};
    end
    Burstloc.bursttime_fit_table = array2table(Burstloc.burst_time_fit, 'VariableNames', Burstloc.col_var_name, 'RowNames', Burstloc.row_var_name);
    writetable(Burstloc.bursttime_fit_table,'Burst locomotion time fit.csv', 'WriteRowNames', true, 'Delimiter', 'tab');
    Burstloc.bursttime_w_table = array2table(Burstloc.burst_time_w, 'VariableNames', Burstloc.col_var_name, 'RowNames', Burstloc.row_var_name);
    writetable(Burstloc.bursttime_w_table,'Burst locomotion time window.csv', 'WriteRowNames', true, 'Delimiter', 'tab');
    
    figure; hold on
    subplot(2,1,1);
    x = 1:size(Burstloc.bursttime_fit_table,2);
    y = Burstloc.bursttime_fit_table{:,:};
    plot(x, y, '.', 'MarkerEdgeColor', 'b', 'MarkerSize', 30);
    title('onset fit detection');
    xlabel('Axon');
    ylabel('Onset (ms)');
    subplot(2,1,2);
    x = 1:size(Burstloc.bursttime_w_table,2);
    y = Burstloc.bursttime_w_table{:,:};
    plot(x, y, '.', 'MarkerEdgeColor', 'b','MarkerSize', 30);
    title('onset window detection');
    xlabel('Axon');
    ylabel('Onset (ms)');
    hold off
    saveas(gcf, 'burst locomotion onset.fig');

end


%% Change of dFF during puff experiment organized in trials.

if Flag_puff == 1;
    x_puff = {};
    y_puff = {};
    TimeAxon_trial_sec = TimeAxon_trial; % convert the time axis in second for the figures.
    for i = 1:Numb_trials;
        TimeAxon_trial_sec{1,i} = TimeAxon_trial_sec{1,i}./1e3;
    end
    prompt = {'air puff duration (in sec)','baseline (in sec)'};
    dlg_title = 'air puff';
    num_lines = 1;
    defaultans = {'','2'};
    ans_puff = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Puff_duration = str2num(ans_puff{1,1});
    Baseline_puff = str2num(ans_puff{2,1});
    
    % figure with traces for each axon and puff highlighted.
    y_offset = 1;
    figure; hold on
    for i = 1:size(Axon_dFF_smooth, 2);
        Y_offset = (y_offset * i) - y_offset;
        plot(TimeAxon(:,1)./1000, Axon_dFF_smooth(:,i) + Y_offset);
    end
    Trial_length = (TimeAxon(end,1)/1000)/Numb_trials; % in sec.
    xlabel('Time (s)');
    axis([0 TimeAxon(end,1)/1000 -2 inf]);
    hold off
    hold on
    for i = 1:Numb_trials;
        rectangle('position',[Baseline_puff+(Trial_length*i-Trial_length) -2 Puff_duration y_offset*Numb_total_axon+2], 'EdgeColor', 'g', 'LineWidth', 1);
    end
    title('dFF smooth all axons');
    saveas(gcf, 'dFF smooth all axons and puffs.fig');
    hold off
    
    dFF_smooth_trials_norm = Axon_dFF_smooth_trial;
    for i = 1:Numb_trials;
        figure; hold on;
        dFF_smooth_trials_norm{1,i} = bsxfun(@rdivide, Axon_dFF_smooth_trial{1,i}, max(Axon_dFF_smooth_trial{1,i}));
        subplot(3,1,1);
        imagesc(dFF_smooth_trials_norm{1,i}', [0 1]);
        colormap(jet);
        axis([0 size(dFF_smooth_trials_norm{1,i}, 1) 0 size(dFF_smooth_trials_norm{1,i}, 2)+1]);
        axis ij
        subplot(3,1,2);
        x_puff{1,i} = TimeAxon_trial_sec{1,i}(:,1)-TimeAxon_trial_sec{1,i}(1,1);
        y_puff{1,i} = zeros(size(TimeAxon_trial_sec{1,i},1),1);
        [row,~] = find((Baseline_puff < x_puff{1,i}(:,1) & x_puff{1,i}(:,1) < Baseline_puff+Puff_duration));
        y_puff{1,i}(row,1) = 0.5;
        plot(x_puff{1,i}, y_puff{1,i}, 'g', 'LineWidth', 1);
        axis([x_puff{1,i}(1,1) x_puff{1,i}(end,1) 0 1]);
        %set(gca,'visible', 'off');
        subplot(3,1,3);
        if Loco_data == 1;
            plot(SpeedTimeTrial_sec{1,i}, SpeedDataTrial{1,i}, 'k', 'LineWidth', 1);
            axis([SpeedTimeTrial_sec{1,i}(1,1) SpeedTimeTrial_sec{1,i}(end,1) -5 inf]);
            xlabel('time (s)');
        end
        saveas(gcf, sprintf('PFs puff and locomotion trial_%d.fig',i));
        hold off
    end
    
    % plot dFF before, during and after the end of the puff.
    
    for i = 1:Numb_trials;
        figure; hold on
        y_offset = 0.1;
        subplot(2,1,1);
        for ii = 1:Numb_total_axon;
            hold on
            Y_offset = (y_offset * ii) - y_offset;
            [row, ~] = find(TimeAxon_trial{1,i}(:,ii)./1000 > ((Baseline_puff+TimeAxon_trial{1,i}(1,ii)/1000)-0.5) & TimeAxon_trial{1,i}(:,ii)./1000 < (TimeAxon_trial{1,i}(1,ii)/1000)+(Baseline_puff+Puff_duration+3)); % values are adjusted to be in seconds.
            x = 1:size(row,1);
            plot(x', Axon_dFF_smooth_trial{1,i}(row,ii) + Y_offset);
            hold off
        end
        axis([0 inf -0.5 inf]);
        subplot(2,1,2);
        [row, ~] = find(x_puff{1,i}(:,1) > Baseline_puff-0.5 & x_puff{1,i}(:,1) < Baseline_puff+Puff_duration+3);
        plot(x_puff{1,i}(row,1), y_puff{1,i}(row,1), 'g', 'LineWidth', 1); % time in sec.
        axis([x_puff{1,i}(row(1,1),1) x_puff{1,i}(row(end,1),1) 0 1]);
        saveas(gcf, sprintf('PFs puff zoom in trial %d.fig',i))
        hold off
    end
    
    % quantify timing of bursts at the beginning and end of the puff.
    a = exist('Frame_rate');
    if a == 1;
        acquisition_rate = Frame_rate;
    end
    clear a
    clear Row_max K_w X_f fit_function burst_time_threshold burst_time_fit burst_time_w fit_function
    Burstpuff.Burst_time_threshold = zeros(Numb_trials,Numb_total_axon); % create a matrix that will contain the time of the first burst when it reaches the threshold.
    Burstpuff.burst_time_fit = NaN(Numb_trials,Numb_total_axon); % create a matrix that will contain the time of the first burst during the puff for each axons.
    Burstpuff.burst_time_w = zeros(Numb_trials,Numb_total_axon);
    Burstpuff.row_max = NaN(Numb_trials,Numb_total_axon);
    Burstpuff.row_event = zeros(Numb_trials,Numb_total_axon);
    Burstpuff.fit_function = cell(Numb_trials,Numb_total_axon);
    Burstpuff.onset = zeros(Numb_trials,Numb_total_axon);
    Burstpuff.X_f = cell(Numb_trials,Numb_total_axon);
    Burstpuff.K_w = cell(Numb_trials,Numb_total_axon);
    
    
    for i = 1:Numb_total_axon;
        for ii = 1:Numb_trials;
            %[burst_time_threshold, Row_event, Row_max, Fit_function, Onset, X, K,...
            %burst_onset_f, burst_onset_w] = burst_detection_trigger(data,...
            %data_time, acquisition_rate, trigger_time, baseline, detection_window);
            %values are in ms.
            [burst_time_threshold, Row_event, Row_max, Fit_function, Onset, X, K, burst_onset_f, burst_onset_w] = burst_detection_trigger(Axon_dFF_smooth_trial{1,ii}(:,i), TimeAxon_trial{1,ii}(:,i), acquisition_rate, Baseline_puff*1000, 500, 3000);
            Burstpuff.Burst_time_threshold(ii,i) = burst_time_threshold;
            Burstpuff.burst_time_fit(ii,i) = burst_onset_f;
            Burstpuff.burst_time_w(ii,i) = burst_onset_w;
            Burstpuff.row_max(ii,i) = Row_max;
            Burstpuff.row_event(ii,i) = Row_event;
            Burstpuff.fit_function{ii,i} = Fit_function;
            Burstpuff.onset(ii,i) = Onset;
            Burstpuff.X_f{ii,i} = X;
            Burstpuff.K_w{ii,i} = K;
        end
    end
    
    % figures showing where the onset has been taken.
    prompt = {'check traces for onset?'};
    dlg_title = 'yes=1, no=0';
    num_lines = 1;
    defaultans = {'1'};
    ans_onset = inputdlg(prompt,dlg_title,num_lines,defaultans);
    ans_onset = str2num(ans_onset{1,1});
    if ans_onset == 1;
        close all
        for i = 1:Numb_total_axon;
            for ii = 1:Numb_trials;
                flag_w = 0;
                flag_fit = 1;
                a = isfield(Burstpuff, 'fit_function');
                if a == 1;
                    f = Burstpuff.fit_function{ii,i};
                    flag_fit = isempty(f);
                    x = Burstpuff.X_f{ii,i};
                end
                b = isfield(Burstpuff, 'K_w');
                if b == 1
                    k = Burstpuff.K_w{ii,i};
                    if k > 0;
                        flag_w = 1;
                    end
                end
                subplot(2,1,1);
                if flag_fit == 0;
                    plot(f, x', Axon_dFF_smooth_trial{1,ii}(Burstpuff.row_event(ii,i):Burstpuff.row_max(ii,i),i), 'k');
                    if round(Burstpuff.onset(ii, i)) > 0 & round(Burstpuff.onset(ii, i)) < size(Axon_dFF_smooth_trial{1,ii},1);
                        hold on
                        plot(round(Burstpuff.onset(ii,i)), Axon_dFF_smooth_trial{1,ii}(Burstpuff.row_event(ii,i)+round(Burstpuff.onset(ii,i)),i), 'r*');
                        hold off
                    end
                end
                title(['fit detection', sprintf('Axon %d', i), sprintf('Trial %d',ii)]);
                subplot(2,1,2);
                if row_max(ii,i) > 0;
                    plot(Axon_dFF_smooth_trial{1,ii}(Burstpuff.row_event(ii,i):Burstpuff.row_max(ii,i),i), 'k');
                    if flag_w == 1;
                        hold on
                        plot(k-Burstpuff.row_event(ii,i)+1, Axon_dFF_smooth_trial{1,ii}(k,i), 'r*');
                        hold off
                    end
                end
                title('sliding window detection');
                clear f x k flag_fit a b
                waitfor(gcf);
            end
        end
    end
    clear flag_w
    
    % burst time is baseline substracted.
    Burstpuff.burst_time_fit = Burstpuff.burst_time_fit - Baseline_puff*1000; % in ms, baseline substracted.
    for i = 1:(size(Burstpuff.burst_time_fit,1)); % remove operation performed on NaN values.
        for ii = 1:(size(Burstpuff.burst_time_fit,2));
            if Burstpuff.burst_time_fit(i,ii) == -Baseline_puff*1000
                Burstpuff.burst_time_fit(i,ii) = NaN;
            end
        end
    end
    Burstpuff.burst_time_w = Burstpuff.burst_time_w - Baseline_puff*1000; % in ms, baseline substracted.
    for i = 1:(size(Burstpuff.burst_time_w,1)); % remove operation performed on NaN values.
        for ii = 1:(size(Burstpuff.burst_time_w,2));
            if Burstpuff.burst_time_w(i,ii) == -Baseline_puff*1000
                Burstpuff.burst_time_w(i,ii) = NaN;
            end
        end
    end
    % creates table containing the burst data.
    for i = 1:size(Burstpuff.burst_time_fit,2);
        Burstpuff.col_var_name(1,i) = {sprintf('axon_%d',i)};
    end
    for i = 1:Numb_trials;
        Burstpuff.row_var_name(1,i) = {sprintf('trial_%d',i)};
    end
    Burstpuff.bursttime_fit_table = array2table(Burstpuff.burst_time_fit, 'VariableNames', Burstpuff.col_var_name, 'RowNames', Burstpuff.row_var_name);
    writetable(Burstpuff.bursttime_fit_table,'Burst puff time fit.csv', 'WriteRowNames', true, 'Delimiter', 'tab');
    Burstpuff.bursttime_w_table = array2table(Burstpuff.burst_time_w, 'VariableNames', Burstpuff.col_var_name, 'RowNames', Burstpuff.row_var_name);
    writetable(Burstpuff.bursttime_w_table,'Burst puff time window.csv', 'WriteRowNames', true, 'Delimiter', 'tab');
    
    figure; hold on
    subplot(2,1,1);
    x = 1:size(Burstpuff.bursttime_fit_table,2);
    y = Burstpuff.bursttime_fit_table{:,:};
    plot(x, y, '.', 'MarkerEdgeColor', 'b', 'MarkerSize', 30);
    title('onset fit detection');
    xlabel('Axon');
    ylabel('Onset (ms)');
    subplot(2,1,2);
    x = 1:size(Burstpuff.bursttime_w_table,2);
    y = Burstpuff.bursttime_w_table{:,:};
    plot(x, y, '.', 'MarkerEdgeColor', 'b','MarkerSize', 30);
    title('onset window detection');
    xlabel('Axon');
    ylabel('Onset (ms)');
    hold off
    saveas(gcf, 'burst puff onset.fig');
    
end


%% events detection. By using event_detection2 function.
% function [detection_criterion, peaks, loc, fitted_template] =...
% event_detection2(data, duration, rise_time, decay_time, threshold, make_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% needs to be updated from what
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% has been done in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROI_grouping_2.mat
if Flag_eventdetection == 1;
    
    clear Std p Loc Peaks
    
    Rise_time = 50; % in ms.
    Decay_time = 100; % in ms.
    p = {}; % dFF max values for the peak detected.
    Std = {}; % standard deviation used for threshold peak detection.
    Loc = {}; % x values of the peaks.
    
    for i = 1:Numb_total_axon; % perform the analysis on each axon.
        
        SD = std(Axon_z_smooth(10000:14000,i)); % SD calculated on a portion of trace where there is not much activity.
        [~, peaks, loc, ~] = event_detection2(Axon_z_smooth(:,i)', TimeTotal_exp, Rise_time, Decay_time, 5 * SD, false);
        % threshold set at 'x' times the SD.
        Std{i} = SD;
        p{i} = peaks;
        Loc{i} = loc;
        
    end
    
    Peaks = {}; % give the dFF peak values from the matrix data.
    for i = 1:Numb_total_axon; % perform the analysis on each axon.
        % plot the detected peak on the dFF traces and make a figure for each
        % axon.
        Peaks{i} = Axon_dFF_smooth(Loc{1,i}(1,:),i);
        
        figure; % make figures for all the calcium traces.
        plot(Axon_dFF_smooth(:,i));
        hold on
        plot(Loc{1,i}(1,:), Peaks{1,i}(:,1), 'ro');
        hold off
        
    end
end

