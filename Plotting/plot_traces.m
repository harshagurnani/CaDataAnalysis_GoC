function [] = plot_traces( data, times, tp_roi_trial)
%% Syntax : plot_traces( data, times, [tp, roi, trial] )
% Given: 'data' with dim index array 'tp_roi_trial' ([tp, roi, trial]), plot traces : POIs on diff rows/ trials one after the other along the row
% ROI is the single dimension of space - could be POI, averaged patches or
% volumes (NN --> N)
% ftype - Scaling ptions:
%         'dFF' --> between 0 and 1
%         'roi' --> individual poi has its own scale - between its own min/max

    %% Sixe of data
    cmap = colormap(lines);
    tp = size(data,tp_roi_trial(1));
    
    %How many trials
    if ~isnan(tp_roi_trial(3))
        trial = size(data,tp_roi_trial(3));
    else
        trial=1;
    end
    
    
    %How many ROIs (for per ROI normalization)
    if ~isnan(tp_roi_trial(2))
        roi = size(data,tp_roi_trial(2));
    else
        roi=1;
    end
    
    for tn = 1:trial
        for rn = 1:roi
            plot(times(:,tn,rn), data(:,tn,rn)+ roi-rn+1, 'Color',cmap(mod(tn,7)+1,:))
            hold on
        end
    end
end

