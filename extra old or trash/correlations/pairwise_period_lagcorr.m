function [lags, corrs, whisk_period] =  pairwise_period_lagcorr( POI, dff, dff_time, ds_rate, lag_halfrange, varargin )
%% SYNTAX: [l, c] = pairwise_period_lagcorr( locoperiod, dff_values, dff_times,  Downsampling_rate, MaxLag [, extratime )
%   Computes pairwise Correlation coefficient between dff for all ROIs
%   Correlation is computed after both data has been cropped to respective
%   periods, downsampled to 'ds_rate' and within a range of
%   [-lag_halfrange, lag_halfrange] (in ms).
% 
%   ARGUMENTS:
%       - POI   [2-D array]
%           Periods Of Interest. 
%           Start and end times (ms) of periods of interest should be first
%           two columns.
%       - dff    [2-D array]
%           All dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - dff_time    [2-D array]
%           All Timestamps (ms) for dFF data. Dim 1 is time, Dim 2 is ROI number.
%       - ds_rate   [INT]
%           Downsampling rate for computing correlations (in Hz)
%       - lag_halfrange     [INT]
%           Maximum lag(or lead) (in ms)  for cropping correlations
%
%   (optional)
%       - Varargin{1}: extratime     [INT or [INT,INT]]
%           Extratime (in ms) around period of interest.
%   OUTPUT:
%       - lags [1-D array]
%           Lag values (in ms)
%       - corrs [4-D array]
%           Lag correlation values ('coeff'). Dim 1 & 2 is ROI number, Dim
%           3 is lag value (in ms). Dim 4 is period of interest.
%       - pkperiod_corr [ 3-D array]
%           Peak pairwise corrn for each period. Dim 1 & 2 is ROI number,
%           Dim 3 is period of interest.
%       - pkcorr [ 2-D array]
%           Peak pairwise corrn. Dim 1 & 2 is ROI number.
%
%
%%%%%%%%% TO ADD!!!!
%%%%%%%%% Allow inputs to be cell arrays instead of vectors
%%%%%%%%% Allow sending cropped data
%%%%%%%%% Different modes of correlations

    %% Parse inputs
    extratime = []; min_POI_length = 0; do_plot = true;
    if nargin > 7, extratime = varargin{1};         end
    if nargin > 8, min_POI_length = varargin{2};    end
    if nargin > 9, do_plot = varargin{3};           end
    
    if isempty(extratime),      extratime = 0;      end
    if isempty(min_POI_length), min_POI_length = 0; end

    %Corr only between +/- 1.5 s BY DEFAULT
    if isempty(lag_halfrange), lag_halfrange = 1500;    end
    
    whisk_period = [];
    %% Which whisking periods to analyse?
    for kk=1:size(POI,1)
        if POI(kk,2) - POI(kk,1) >= min_POI_length  %Minimum size of POI
            whisk_period=[whisk_period; kk];
        end
    end

    n_periods = max(size(whisk_period));    %No. of whisking periods
    n_ROIs = size(dff,2);   %No. of ROIs

    lags  =(-lag_halfrange:1000/ds_rate:lag_halfrange);     %ms
    if mod(size(lags,2), 2)==0
        disp('Making number of lags odd, and ensuring lag 0 included')
        lags = (-lag_halfrange-500/ds_rate:1000/ds_rate:lag_halfrange+500/ds_rate);
    end
    n_lags = floor(size(lags,2)/2);

    %Array indexing: corrs( roi1, roi2, lag, whisk_period ) 
    corrs = nan( n_ROIs,  n_ROIs, 2*n_lags+1, n_periods);
    
%     pkperiod_corr = nan( n_ROIs,  n_ROIs, n_periods);
%     pkcorr = nan( n_ROIs,  n_ROIs);
    
    
    tic
    
    for roi =1:n_ROIs
        [ dff_per{roi}, dff_pertime{roi} ] = get_TSOI( dff(:,roi), dff_time(:,roi), POI(whisk_period,:), extratime);
    end
    
    for roi1= 1:n_ROIs
    for roi2 = roi1:n_ROIs 
       for ww=1:max(size(whisk_period))
          try
            [c,l] = get_lagcorr_ds( dff_per{roi1}{ww}, dff_pertime{roi1}{ww}, dff_per{roi2}{ww}, dff_pertime{roi2}{ww}, ds_rate);
            %Mid = zero lag index
            mid=(size(l,2)+1)/2;
          
            %Crop between relevant lags
            if 2*n_lags+1 <= size(c,1)
                corrs(roi1, roi2,:,ww) = c( mid-n_lags:mid+n_lags);
                corrs(roi2, roi1,:,ww) = c( mid-n_lags:mid+n_lags);
            else 
                corrs(roi1, roi2, (n_lags+1)-(mid-1):(n_lags+1)+(mid-1), ww) = c;
                corrs(roi2, roi1, (n_lags+1)-(mid-1):(n_lags+1)+(mid-1), ww) = c;
            end
          catch ME
            disp(sprintf('ROI %d, ROI %d, Period %d \n %s', roi1, roi2, whisk_period(ww), ME.message))
          end
%           indx = find( abs(corrs(roi1, roi2, :, ww)) == max(abs(corrs(roi1, roi2, :, ww))));
%           
%           if ~isempty( indx )
%           pkperiod_corr( roi1, roi2, ww) = corrs(roi1, roi2, indx(1), ww);  %Peak per pair per period
%           pkperiod_corr( roi2, roi1, ww) = corrs(roi1, roi2, indx(1), ww);
%           end
              
       end
%        indx = find( abs(corrs(roi1, roi2, :)) == max(abs(corrs(roi1, roi2, :))));
%        if ~isempty(indx)
%        pkcorr( roi1, roi2) = corrs(roi1, roi2, indx(1));        %Peak per pair
%        pkcorr( roi2, roi1) = corrs(roi1, roi2, indx(1));
%        end
    end
    end
    toc

    
    %% Plotting
    
    if do_plot
        figure()
        colormap jet
        imagesc( sum(corrs(:,:, n_lags+1,:),4)/n_periods, [0,1] )     
        title('Zero-lag correlations')  

        n_rows=3;
        n_cols=3;
        for id=1:n_periods
            if mod(id-1,9) == 0 
                figure();
                colormap jet
                suptitle('Zero-lag correlations in different whisking/locomotion periods')   
            end
            subplot(n_rows,n_cols, mod(id-1,9)+1)
            imagesc( corrs(:,:,n_lags+1, id) , [0,1] )
        end


%         n_rows=3;
%         n_cols=3;
%         for id=1:n_periods
%             if mod(id-1,9) == 0 
%                 figure();
%                 colormap jet
%                 suptitle('Peak correlations in Different periods')  
%             end
%             subplot(n_rows,n_cols, mod(id-1,9)+1 )
%             imagesc( pkperiod_corr(:,:,id), [0,1] )     
%         end
%   
% 
%         figure();
%         colormap jet
%         imagesc( pkcorr(:,:), [0,1] )  
%         title('Peak correlations')  
%     end
    
end

function [ extra, doplot ] = parse_options( nargs, varargs, copts)

extra = 300;
doplot = true;

if nargs > copts  

    charctr=0;
    numctr=0;
    start_read_var = copts+1;
    while start_read_var <= nargs
        vn = start_read_var-copts;
       if ischar( varargs{vn} )
           charctr=1;
           if strcmp(varargs{vn}, 'extratime' )
               extra = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'do_plot' )
               doplot = varargs{vn+1};
               start_read_var = start_read_var + 2;     
           else 
                error('Argument name %s is not an optional argument', varargs{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   extra = varargs{vn};
                   start_read_var = start_read_var+1;
               case 2
                   doplot = varargs{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end
end