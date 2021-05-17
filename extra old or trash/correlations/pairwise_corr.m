function [pkcorr, zerocorr, corrs, lags ] = pairwise_corr( dFF, dFF_time, varargin )
%% Pairwise correlations after downsampling for all traces    
%
% OPTIONAL ARGS:
%   - 'ds_rate' -       Downsampling factor
%   - 'maxlag'  -       Max number of time bins for computing lag correlations
%   - 'do_plot' -       Logical - Plot correlations or not. TRUE by default
Tm = size(dFF,1);
nTrials = size(dFF,2);
nROI = size(dFF,3);

[maxlag, DS_rate, doplot ] = parse_options( nargin, varargin, 2);

temp_dff = nan( ceil(Tm/DS_rate), nTrials, nROI);
temp_time = temp_dff;

%% Filter and downsample
for roi = 1:nROI
    for trial=1:nTrials
        temp_dff(:, trial, roi) = decimate( dFF(:,trial,roi), DS_rate);
        temp_time(:, trial, roi) = downsample(dFF_time(:,trial,roi), DS_rate);
    end
end

%% Calculate correlation
corrs = nan(2*maxlag+1, nROI, nROI);
pkcorr = nan(nROI, nROI);
zerocorr = pkcorr;
for roi1 = 1:nROI
    for roi2 = 1:nROI
       [corrs(:,roi1, roi2), lags] = xcorr( reshape(temp_dff(:,:,roi1),ceil(Tm/DS_rate)*nTrials,1) - mean(reshape(temp_dff(:,:,roi1),ceil(Tm/DS_rate)*nTrials,1)), reshape(temp_dff(:,:,roi2),ceil(Tm/DS_rate)*nTrials,1) - mean(reshape(temp_dff(:,:,roi2),ceil(Tm/DS_rate)*nTrials,1)), maxlag, 'coeff');
       pk = max( abs(corrs(:,roi1,roi2)));
       pks = corrs( abs(corrs(:,roi1,roi2))==pk, roi1, roi2 );
       pkcorr(roi1, roi2) = median(pks);
       zerocorr(roi1,roi2) = corrs(maxlag+1, roi1,roi2);
    end
end

if doplot
    figure()
    imagesc(pkcorr, [-1,1])
    title('Peak correlations')

    figure()
    imagesc(zerocorr, [-1,1])
    title('Zero lag correlations')
    caxis([0 1])
end

end


%-------------------------------------------------------------------------%
%% Parsing options
%-------------------------------------------------------------------------%

function [ maxlag, DS_rate, doplot ] = parse_options( nargs, varargs, copts)
    DS_rate = 5;
    maxlag = 5;     %bins
    doplot = true;
    
    if nargs > copts  

        charctr=0;
        numctr=0;
        start_read_var = copts+1;
        while start_read_var <= nargs
            vn = start_read_var-copts;
           if ischar( varargs{vn} )
               charctr=1;
               if strcmp(varargs{vn}, 'ds_rate' )
                   DS_rate = varargs{vn+1};
                   start_read_var = start_read_var + 2;
               elseif strcmp(varargs{vn}, 'maxlag' )
                   maxlag = varargs{vn+1};
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
                       DS_rate = varargs{vn};
                       start_read_var = start_read_var+1;
                   case 2
                       maxlag = varargs{vn};
                       start_read_var = start_read_var+1;
                   case 3
                       doplot = varargs{vn};
                       start_read_var = start_read_var+1;
                   
               end
           end
        end

    end


end
