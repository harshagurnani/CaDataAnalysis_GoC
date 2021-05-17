function [loc, detection_criterion, peaks, fitted_template, scale] = template_event_detection( data, dt, varargin )
%% Event detection using a scaled template-matching scheme (Clements & Bekkers 1997)
% Find events similar to an exponential filter described with rise and
% decay times. Algorithm fits a scaled-and-offseted version of the template
% to data, and detects events where the SSE of the fit is less
% than threshold for detection criterion. 
% Further work needed on choosing optimal rise and decay times
%
%   Usage:
%
%
%   Inputs:
%
%
%
%
%
%   Outputs:
%
%
%
% Adapted from Geoffrey Evans' code by Harsha G.
% April 2018
    
    [p.tmp, p.events, p.do_plot ] = parse_args(varargin{:});
    if isempty(dt), dt = 1; end

    Tm      = size(data,1);
    nROI    = size(data,2);
    if isscalar(p.events.min_scale), p.events.min_scale = repmat(p.events.min_scale,nROI,1);   end  
    if isempty(p.events.min_scale),  p.events.min_scale = repmat(-Inf,nROI,1);                 end
    
    detection_criterion = zeros(Tm, nROI);
    
    template = get_template(Tm, dt, p.tmp.rise_time, p.tmp.decay_time);   %entire template - same length as data
    if p.events.align_peak, [~,peak_idx] = findpeaks(template); end
    
    scale = nan(Tm, nROI);
    
    switch p.tmp.method
        case 'whole'    % default
            fitted_template = zeros( Tm, Tm, nROI); % template-window_offset-roi
            for roi = 1:nROI
                template = get_template(Tm, dt, p.tmp.rise_time, p.tmp.decay_time);   %entire template - same length as data
            for time = 1:Tm
                template = [NaN, template(1:end-1)];
                scale(time,roi) = get_scale(template', data(:,roi), Tm);
                offset = (nansum(data(:,roi)) - scale(time,roi) .* nansum(template))/Tm;   %only take data in cropped 

                fitted_template(:,time,roi) = template .* scale(time,roi) + offset;
                sse = nansum((data(:,roi) - fitted_template(:,time, roi)).^2);
                standard_error = sqrt(sse / (Tm-1)); 

                detection_criterion(time,roi) = scale(time,roi) ./ standard_error;
            end
                [peaks{roi}, loc{roi}] = get_events(p.events.threshold, detection_criterion(:,roi));
                loc{roi} = loc{roi}(scale(loc{roi})>p.events.min_scale(roi));
                peaks{roi} = detection_criterion(loc{roi},roi);
                if p.events.align_peak, loc{roi} = loc{roi}+peak_idx; end
            end

        case 'cropped'  % uses only a small fixed template and moves it along the window
            
            fitted_template = zeros( p.tmp.window, Tm, nROI); % template-window_offset-roi
            for roi = 1:nROI
                template = get_template(p.tmp.window, dt, p.tmp.rise_time, p.tmp.decay_time);   %entire template - same length as data
                data_temp = [data(:,roi); nan(p.tmp.window,1)];                
            for time = 1:Tm-p.tmp.window+1
                cropped = data_temp(time:time+p.tmp.window-1);
                scale(time,roi) = get_scale(template', cropped, p.tmp.window);
                offset = (nansum(cropped) - scale(time,roi) .* nansum(template))/p.tmp.window;   %only take data in cropped 

                fitted_template(:,time,roi) = template .* scale(time,roi) + offset;
                sse = nansum((cropped - fitted_template(:,time, roi)).^2);
                standard_error = sqrt(sse / (p.tmp.window-1)); 

                detection_criterion(time,roi) = scale(time,roi) ./ standard_error;
            end
                [peaks{roi}, loc{roi}] = get_events(p.events.threshold, detection_criterion(:,roi));
                loc{roi} = loc{roi}(scale(loc{roi},roi)>p.events.min_scale(roi));
                peaks{roi} = detection_criterion(loc{roi},roi);
                if p.events.align_peak, loc{roi} = loc{roi}+peak_idx; end
            end            
            
        case 'weighted'
            
            
            
            
    end
    
    if p.do_plot
        figure;
        hold on
        for roi=1:nROI
        plot((0:Tm-1), detection_criterion(:,roi)+5*roi, 'b')
        plot(loc{roi}, peaks{roi}+5*roi, 'ro')
        end
        hold off
    end
end

function [peaks, loc] = get_events(threshold, criterion)
criterion = [0; criterion];
[peaks, loc] = findpeaks(criterion,'MINPEAKHEIGHT',threshold);
loc=loc-1;
end

function scale = get_scale(template, data, N )
    scale_numerator = nansum(template.*data) - nansum(template).*nansum(data)/N;
    scale_denominator = nansum(template.^2) - nansum(template).^2/N;
    scale = scale_numerator ./ scale_denominator;
end

function template = get_template(N, dt, rise_time, decay_time)
    t = (0:N-1) * dt;
    template = (1 - exp(-t/rise_time)) .* exp(-t/decay_time);
    template = template/max(template);
end

function [template, events, do_plot ] = parse_args(varargin)
    % Using Chen et al 2013, tau_peak = 45 ms, tau_half_decay = 142 ms,
    % calculated the correct parameter values
    % I hate them for making me calculate this 100000000 times
    template.rise_time = 30.5; % ms
    template.decay_time = 102.8;   % ms
    template.method = 'whole';
    template.window = 50;
    template.weight = 1;
    template.given  = [];
    events.threshold = -4;
    events.align_peak = false;
    events.min_scale = [];
    do_plot = false;
    
    v1 = 1;
    while v1<=nargin
        if strcmp( varargin{v1}, 'rise' )
            template.rise_time = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'decay' )
            template.decay_time = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'threshold' )
            events.threshold = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'do_plot' )
            do_plot = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'align_peak' )
            events.align_peak = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'min_scale' )
            events.min_scale = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'temp_method' )
            template.method = varargin{v1+1};
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'temp_window' )
            template.window = varargin{v1+1};
            template.method = 'cropped';
            v1=v1+2;
        elseif strcmp( varargin{v1}, 'use_temp' )
            template.given = varargin{v1+1};
            v1=v1+2;
        end
            
    
    end

    
end