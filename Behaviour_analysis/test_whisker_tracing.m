function [varargout] = test_whisker_tracing( video, whiskers, measurements, frame_to_check, whisker_to_check )
%% Preview traced whiskers in video frames
% draw whiskers on frames to check

% video -           Path to video file#
% whiskers -        Path to .whiskers file/ whiskers struct
% measurements -    Path to .measurements file/ measurements struct
% frame_to_check -  Preview specific frames (or 0 for random selection, -1
%                   for all frames )

    
    n_frames = length( imfinfo(video) );    %number of video frames
    if frame_to_check == 0
        % random selection of frames
        frame_to_check = unidrnd( n_frames, [ ceil(0.005*n_frames), 1] );
    elseif frame_to_check == -1 
        % all frames
        frame_to_check = 1:n_frames;
    end
    nt = length( frame_to_check ); frame_to_check = reshape( frame_to_check, [1, nt]);
    if ~isstruct( whiskers ), whiskers = LoadWhiskers( whiskers ); end
    if ~isstruct( measurements ), measurements = LoadMeasurements( measurements ); end
    
    % sort whiskers traces by id
    n_whisk = 50;
    [whisk_trace, whisk_pos, whisk_time] = deal( cell( n_whisk, 1) );
    for jj=1:n_whisk, whisk_trace{jj}={}; whisk_pos{jj} = []; whisk_time{jj}=[]; end
    whisk_trace_temp = whisk_trace;
    

    for t=1:length(whiskers)
       %sort all
       wid = whiskers(t).id+1;
       frame = whiskers(t).time+1;
       whisk_trace_temp{wid}{frame} = [whiskers(t).x, whiskers(t).y];
    end

    for t=1:length(measurements)
       % re-classify by whisker label
       whisk_label = measurements(t).label+1;
       wid = measurements(t).wid+1;
       frame = measurements(t).fid+1;
       if whisk_label ~= 0
        whisk_trace{whisk_label}{frame} = whisk_trace_temp{wid}{frame};
        whisk_pos{whisk_label} = [ whisk_pos{whisk_label}; measurements(t).angle, measurements(t).tip_x, measurements(t).tip_y, measurements(t).follicle_x, measurements(t).follicle_y];
        whisk_time{whisk_label} = [whisk_time{whisk_label}; frame];
       end
    end
    clear whisk_trace_temp
    
    % find whiskers to keep
    ff = arrayfun( @(w) length(whisk_time{w}), 1:n_whisk );
    tokeep = ff>5;
    whisk_trace = whisk_trace(tokeep);
    whisk_pos = whisk_pos(tokeep);
    whisk_time = whisk_time(tokeep);
    n_whisk = length(whisk_trace);
    if isempty( whisker_to_check) , whisker_to_check = 1:n_whisk; end
    %% Now plot sequentially
    figure;
    cmap = colormap( jet(n_whisk) );
    colormap(gray);
    for t=frame_to_check
        % Plotting 8-bit images
        imagesc( imread( video, t ),[150 220] );
        hold on;
        % Plot all whiskers
        for wid=1:n_whisk
           if ~isempty( find(whisker_to_check==wid,1))
           tid = find( whisk_time{wid}==t,1 );
           if ~isempty(tid)
               plot( whisk_trace{wid}{t}(:,1), whisk_trace{wid}{t}(:,2), 'color', cmap(wid,:) );
%                scatter( whisk_pos{wid}(tid, 2), whisk_pos{wid}(tid, 3), 'filled', 'markerfacecolor', cmap(wid+1,:) );
           end
           end
        end
        hold off;
        waitforbuttonpress;
    end
   
    if nargout>0,   varargout{1} = whisk_trace; end
    if nargout>1,   varargout{2} = whisk_pos;   end
    if nargout>2,   varargout{3} = whisk_time;  end 
    
end