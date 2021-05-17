function [burst_time_threshold, Row_event, Row_max, Fit_function, Onset, X, K, burst_onset_f, burst_onset_w] = burst_detection_trigger(data, data_time, acquisition_rate, trigger_time, baseline, detection_window)
% time values should be in ms.
% data are the dFF values and data_time are the time stamp for each point
% in data.
% acquisition_rate in Hz.

Data = data;
Data_time = data_time;

[row_baseline, ~] = find(Data_time(:,1) < trigger_time+Data_time(1,1) & Data_time(:,1) > trigger_time+Data_time(1,1)-baseline); % baseline window before the trigger.
SD_baseline = std(Data(row_baseline,1)); % SD for the baseline.
median_baseline = median(Data(row_baseline,1)); % median for the baseline.
[row_event, ~] = find(Data_time(:,1) > trigger_time+Data_time(1,1) & Data_time(:,1) < trigger_time+Data_time(1,1)+detection_window); % define the window for the event detection.
[row_time, ~] = find(Data(row_event,1) > median_baseline+5*SD_baseline); % threshold set at 5 times the SD before the trigger above the baseline.
Row_event = row_event(1,1);
TF = isempty(row_time);
if TF == 1;
    burst_time_threshold = NaN; % if no burst detected in the window.
    Fit_function = []; % fit_function is empty.
    Onset = NaN;
    X = NaN; % x values for fit_function empty.
    K = NaN; % row value for window detection empty.
    Row_max = NaN; % row value for fit_function empty.
    burst_onset_f = NaN;
    burst_onset_w = NaN;
end
if TF == 0;
    burst_time_threshold = Data_time(row_event(1,1)+row_time(1,1),1)-Data_time(1,1); % find the time in ms of the first burst from the trigger time.
    % burst detection by using a fit.
    for j = row_time(1,1):row_event(1,1)+row_time(1,1)+acquisition_rate % define the window where the max of the burst is found with a 1 second window.
        if Data(row_event(1,1)+j,1) > Data(row_event(1,1)+j+10,1); % 10 points sliding window.
            [~, M] = max(Data(row_event(1,1)+j:row_event(1,1)+j+10,1));
            row_max = row_event(1,1)+j+M;
            break
        end
    end
    Row_max = row_max;
    x = 1:size(Data(row_event(1,1):row_max),1); % create x for the fit function.
    f = fit(x',Data(row_event(1,1):row_max,1)-(2*SD_baseline+median_baseline),'poly3'); % fit the trace from beginning of event window to the detected threshold. The threshold has been substracted from the fit which enable to find the x value for threshold set at 2 times the SD.
    X = x;
    Fit_function = f;
    onset = fzero(f,f.p4); % row values where the event starts in the event detection window. f.p4 is the result for y=0.
    Onset = onset;
    if onset < 0;
        burst_onset_f = NaN;
    end
    if onset > 0 && onset < size(Data_time,1);
        burst_onset_f = Data_time(row_event(1,1)+round(onset),1)-Data_time(1,1); % time where the events start in ms.
    else
        burst_onset_f = NaN;
    end
    % burst detection by using a backward sliding window.
    baseline_w = round((200 * acquisition_rate)/1000); % number of point for a 200 ms baseline window.
    for k = row_time(1,1)+row_event(1,1):-1:row_event(1,1); % search backward.
        SD_w = std(Data(k-baseline_w:k,1));
        median_w = median(Data(k-baseline_w:k,1));
        if Data(k,1) < 2*SD_w+median_w; % threshold set at 2 times the SD.
            burst_onset_w = Data_time(k,1)-Data_time(1,1); % time where the events start in ms.
            K = k;
            break
        end
        if k == row_event(1,1);
            burst_onset_w = Data_time(row_event(1,1),1)-Data_time(1,1);
            K = row_event(1,1);
        end
    end
end

end