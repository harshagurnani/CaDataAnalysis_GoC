%% Detect whisking periods based on the motion index of the pad (Whiskers.MI_whiskerpad) and the whiskers
% angle if analyzed.

% define period of pad motion.

Whiskers.MI_whiskerpad_smooth = smooth(Whiskers.MI_whiskerpad(:,1),10);
figure; plot(Whiskers.MI_whiskerpad_smooth(:,1));
waitfor(gcf);
prompt = {'threshold'};
dlg_title = 'threshold for MI';
num_lines = 1;
defaultans = {''};
ans_MI = inputdlg(prompt,dlg_title,num_lines,defaultans);

% mean_w = mean(Whiskers.MI_whiskerpad(str2num(ans_MI{1,1}):str2num(ans_MI{2,1}),1));
% SD_w = std(Whiskers.MI_whiskerpad(str2num(ans_MI{1,1}):str2num(ans_MI{2,1}),1));

[Whiskers.rowMI, ~] = find(Whiskers.MI_whiskerpad(:,1) > str2num(ans_MI{1,1})); % find values.

if isempty(Whiskers.rowMI) == 0;
    Whiskers.MI_time = Whiskers.MI_whiskerpad(Whiskers.rowMI, 2); % time in ms.
    counter = 1;
    Whiskers.MI_range = {};
    for i = 1:size(Whiskers.MI_time, 1)-1; % find the different whisker pad motion period.
        if Whiskers.MI_time(i+1) > Whiskers.MI_time(i) + 1000; % periods spaced by a minimum of 1000 ms.
            if counter == 1;
                Whiskers.MI_range{counter}(1,1) = Whiskers.MI_time(1);
                Whiskers.MI_range{counter}(1,2) = Whiskers.MI_time(i);
                ii = Whiskers.MI_time(i+1);
            else
                Whiskers.MI_range{counter}(1,1) = ii;
                Whiskers.MI_range{counter}(1,2) = Whiskers.MI_time(i);
                ii = Whiskers.MI_time(i+1);
            end
            counter = counter+1;
        end
    end
    if counter > 1; % find the last motion period.
        Whiskers.MI_range{counter}(1,1) = ii;
        Whiskers.MI_range{counter}(1,2) = Whiskers.MI_time(end,1);
    end
    if isfield(Whiskers, 'MI_range') == 0; % if only one period.
        Whiskers.MI_range{counter}(1,1) = Whiskers.MI_time(1,1);
        Whiskers.MI_range{counter}(1,2) = Whiskers.MI_time(end,1);
    end
    
    Whiskers.Numb_MI_period = size(Whiskers.MI_range, 2);
    
    Whiskers.MI_dur = {};
    for i = 1:size(Whiskers.MI_range, 2);
        Whiskers.MI_dur{i} = (Whiskers.MI_range{1,i}(1,2) - Whiskers.MI_range{1,i}(1,1)) / 1000; % in s.
    end
    Whiskers.MI_duration = zeros(size(Whiskers.MI_range));
    for i = 1:size(Whiskers.MI_range, 2);
        Whiskers.MI_duration(1,i) = Whiskers.MI_dur{1,i}(1,1);
    end
    Whiskers.MI_mean_dur = mean(Whiskers.MI_duration, 2); % in s.
    
else
    Whiskers.Numb_MI_period = 0;
end
clear i counter prompt dlg_title num_lines

% figure with the dFF overtime and the motion index for the pad.

figure; hold on
subplot(2,1,1);
imagesc(Axon_dFF_smooth_norm', [0 1]);
ax = gca;
ax.XAxis.Visible = 'off';
colormap(jet);
if size(Axon_dFF_smooth_norm, 2)>1;
    axis([0 size(Axon_dFF_smooth_norm, 1) 1 size(Axon_dFF_smooth_norm, 2)]);
end
axis ij;
subplot(2,1,2);
plot(Whiskers.MI_whiskerpad(:,2)./1000, Whiskers.MI_whiskerpad(:,1), 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('au');
axis([0 TimeAxon(end,end)./1000 0 inf])
saveas(gcf, 'PFs MI pad.fig');
hold off
clear ax

if Whiskers.Numb_MI_period > 0;
    
    % figure with MI periods and dFF traces.
    y_offset = 0.5;
    figure; hold on
    for i = 1:size(Axon_dFF_smooth, 2);
        Y_offset = (y_offset * i) - y_offset;
        plot(TimeAxon(:,1)./1000, Axon_dFF_smooth(:,i) + Y_offset);
    end
    axis([0 inf -2 inf]);
    for i = 1:Whiskers.Numb_MI_period;
        p = patch([Whiskers.MI_range{1,i}(1,1)/1000 Whiskers.MI_range{1,i}(1,2)/1000 Whiskers.MI_range{1,i}(1,2)/1000 Whiskers.MI_range{1,i}(1,1)/1000], [-2 -2 y_offset*Numb_total_axon+1 y_offset*Numb_total_axon+1], 'b');
        set(p,'FaceAlpha',0.2);
        set(p, 'EdgeColor', 'none');
    end
    %ax = gca;
    %ax.YAxis.Visible = 'off';
    xlabel('Time (s)');
    title('dFF smooth all axons');
    saveas(gcf, 'dFF smooth all axons and MI pad.fig');
    hold off
    
    Whiskers.baseline_MI = 3000; % in ms, time before and after a motion period.
    
    for i = 1:Whiskers.Numb_MI_period;
        for j = 1:Numb_total_axon;
            [row, ~] = find(TimeAxon(:,j) > (Whiskers.MI_range{1,i}(1,1))-Whiskers.baseline_MI & TimeAxon(:,j) < (Whiskers.MI_range{1,i}(1,2))+Whiskers.baseline_MI);
            Whiskers.Row{i,j} = row;
            Axon_dFF_MI{i,j} = Axon_dFF(Whiskers.Row{i,j},j);
            Axon_dFF_smooth_MI{i,j} = Axon_dFF_smooth(Whiskers.Row{i,j},j);
        end
    end
    Whiskers.MI_dim = zeros(size(Axon_dFF_MI)); % all the traces are reshape to have the same length.
    for i = 1:Whiskers.Numb_MI_period;
        for j = 1:Numb_total_axon;
            Whiskers.MI_dim(i,j) = size(Axon_dFF_MI{i,j},1);
        end
    end
    for i = 1:Whiskers.Numb_MI_period;
        M = min(Whiskers.MI_dim(i,:));
        for j = 1:Numb_total_axon;
            Axon_dFF_MI{i,j} = Axon_dFF_MI{i,j}(1:M,1);
            Axon_dFF_smooth_MI{i,j} = Axon_dFF_smooth_MI{i,j}(1:M,1);
        end
    end
    
    % make one figure for each locomotion period. norm for the higher dFF
    % per trace.
    for i = 1:Whiskers.Numb_MI_period;
        figure; hold on
        Whiskers.dFF_smooth_MI = zeros(size(Axon_dFF_smooth_MI{i,1},1),Numb_total_axon);
        for j = 1:Numb_total_axon;
            Whiskers.dFF_smooth_MI(:,j) = Axon_dFF_smooth_MI{i,j}(:,1);
        end
        Whiskers.dFF_smooth_MI_norm = Whiskers.dFF_smooth_MI;
        Whiskers.dFF_smooth_MI_norm = bsxfun(@rdivide, Whiskers.dFF_smooth_MI, max(Axon_dFF_smooth)); % normalized for the max value over the recording.
        subplot(2,1,1);
        imagesc(Whiskers.dFF_smooth_MI_norm', [0 1]);
        colormap(jet);
        ax = gca;
        ax.XAxis.Visible = 'off';
        %axis([0 size(dFF_smooth_loc_norm, 1) 1 size(dFF_smooth_loc_norm, 2)+1])
        axis ij;
        subplot(2,1,2);
        [row, ~] = find(Whiskers.MI_whiskerpad(:,2) > (Whiskers.MI_range{1,i}(1,1))-Whiskers.baseline_MI & Whiskers.MI_whiskerpad(:,2) < (Whiskers.MI_range{1,i}(1,2))+Whiskers.baseline_MI);
        % 2000 ms before and after the locomotion period.
        plot(Whiskers.MI_whiskerpad(row,2)./1000, Whiskers.MI_whiskerpad(row, 1), 'k', 'LineWidth', 1);
        axis([Whiskers.MI_whiskerpad(row(1,1),2)./1000 inf 0 inf]);
        xlabel('Time (s)');
        saveas(gcf, sprintf('MI%d.fig',i));
        hold off
    end
    clear ax
    
    % plot dFF traces during each locomotion period for all axons.
    
    for i = 1:Whiskers.Numb_MI_period;
        figure; hold on
        y_offset = 0.1;
        subplot(2,1,1);
        for ii = 1:Numb_total_axon;
            hold on
            Y_offset = (y_offset * ii) - y_offset;
            x = 1:size(Axon_dFF_smooth_MI{i,ii}(:,1),1);
            plot(x', Axon_dFF_smooth_MI{i,ii}(:,1) + Y_offset);
            axis([0 inf -0.2 inf]);
            hold off
        end
        subplot(2,1,2);
        hold on
        [row, ~] = find(Whiskers.MI_whiskerpad(:,2) > (Whiskers.MI_range{1,i}(1,1))-Whiskers.baseline_MI & Whiskers.MI_whiskerpad(:,2) < (Whiskers.MI_range{1,i}(1,2))+Whiskers.baseline_MI);
        x = (Whiskers.MI_whiskerpad(row, 2)./1000)-(Whiskers.MI_whiskerpad(row(1,1), 2)./1000);
        plot(x', Whiskers.MI_whiskerpad(row, 1), 'k', 'LineWidth', 1);
        axis([0 inf 0 inf]);
        xlabel('Time (s)');
        ylabel('au');
        hold off
        saveas(gcf, sprintf('PFs MI period_%d.fig',i));
        hold off
        clear x row
    end
end



