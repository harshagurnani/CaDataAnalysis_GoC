function [] = plot_scaled_image( data, tp_roi_trial, ftype)
%% Syntax:  plot_scaled_image( data, [tp, roi, trial], ftype )
% Given: 'data' with dim index array 'tp_roi_trial' ([tp, roi, trial]) and scaling type, plot colour-scale map with separate rows for each POI and trials one after the other along the row
% ROI is the single dimension of space - could be POI, averaged patches or
% volumes (NN --> N) ftype - Scaling ptions:
%         'dFF'     --> between 0 and max dFF
%         'roi'     --> individual poi has its own scale - between its own min/max
%         'zscore'  --> Single scale for z-scored data

    %% Reshaping data
    temp = permute(data,tp_roi_trial);
    tp = size(data,tp_roi_trial(1));
    
    Newtemp=[];
    %Concatenate trials
    if ~isnan(tp_roi_trial(3))
        trial = size(data, tp_roi_trial(3));
        for tn = 1:trial
            Newtemp = vertcat(Newtemp,temp(:,:,tn));
            %Newtemp = [Newtemp; zeros(100,roi)];
        end
    else
        trial = 1;
    end
    
    %How many ROIs (for per ROI normalization)
    if ~isnan(tp_roi_trial(2))
        roi = size(data,tp_roi_trial(2));
    else
        roi=1;
    end
    
    clear temp
    Newtemp = transpose(Newtemp);
    
    %% Setting color scale
    switch ftype
        case 'dFF'
            clims = [0 max(Newtemp(:))];
            colormap(parula)
        case 'roi'
            for jj = 1:roi
                m = min(Newtemp(jj,:));
                M = max(Newtemp(jj,:));
                Newtemp(jj,:)= (Newtemp(jj,:) - m) ./(M-m);
                clims = [0 1];
            end
        case 'zscore'
            for jj = 1:roi
               mean_n = mean(Newtemp(jj,:));
               sd_div = std(Newtemp(jj,:)); 
               Newtemp(jj,:) = (Newtemp(jj,:) - mean_n)./sd_div;
               ulim= ceil(max(abs(Newtemp(:))));
               llim= floor(min(Newtemp(:)));
               clims = [-ulim ulim];
            end
            colormap(jet)
    end
    
    imagesc(Newtemp, clims)
end