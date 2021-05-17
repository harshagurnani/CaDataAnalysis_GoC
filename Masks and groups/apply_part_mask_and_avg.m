function [ masked_data, time_mask] = apply_part_mask_and_avg( data, times, cropping, masks, varargin)
%% Apply masks to patches given by roi numbers (so you can use multiple masks in same patch). Normalisation and averaging groups is optional.
% 
% ARGUMENTS:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%

NY = size(data,1);
NX = size(data,2);
Tm = size(data,3);
nTrials = size(data,4);
nr = size(data,5);

%% Parse optional arguments
[ roi_indx, doavg, avg_groups, donorm, params, doplot ] = parse_options( nargin, varargin, 4, nr );
nROI = length(roi_indx) ;

masked_data = nan( Tm, nTrials, nROI);
time_mask = masked_data;

if size(cropping,2) == 1
    DX = cropping;
    DY = cropping;
else
    DX = cropping(2);
    DY = cropping(1);
end

xcrop = DX+1:NX-DX;
ycrop = DY+1:NY-DY;

nx = size(xcrop,2);
ny = size(ycrop,2);



%% Apply masks
for roi=1:nROI
    patch_id = roi_indx(roi);
    temp = reshape(data(ycrop,xcrop,:,:,patch_id), [ny*nx, Tm, nTrials]);
    allx = floor(find( masks{roi} == 1 )/ny)+1;
    ally = mod(find( masks{roi} == 1 ),ny)+1;
    
    c_time = [floor(mean(ally)) floor(mean(allx))];
    
    for trial=1:nTrials
    for t=1:Tm
        masked_data(t,trial,roi) = mean( temp(masks{roi}(:), t, trial), 1 );
        time_mask(t,trial,roi) = times(ycrop(c_time(1)), xcrop(c_time(2)), t, trial, patch_id);
    end
    end
end

%% Optional averaging
if doavg
    ngroups = numel(avg_groups);
    avg_data = nan( Tm, nTrials, ngroups);
    avg_time = avg_data;

    for group=1:ngroups
        % Mean or sum?
        avg_data(:,:,group) = mean( masked_data(:,:, avg_groups{group}),3 );
        avg_time(:,:,group) = mean( time_mask(:,:, avg_groups{group}),3 );
    end
    
    masked_data = avg_data;
    time_mask = avg_time;
    
end
%% Optional normalisation and plotting
if donorm
    [norm, ~] =normalize_and_smooth( masked_data, time_mask, params);
    if doplot
        plot_traces(norm, time_mask, [1 3 2])
    end
    masked_data = norm;
end

end





%-------------------------------------------------------------------------%
%% PARSING OPTIONS %%
%-------------------------------------------------------------------------%

function [ roi_indx, doavg, avg_groups, donorm, params, doplot ] = parse_options( nargs, varargs, copts, nroi)

roi_indx = 1:nroi;
doavg=false;
avg_groups = [];
donorm = [];
params = [];
doplot =[];

if nargs > copts  

    charctr=0;
    numctr=0;
    start_read_var = copts+1;
    while start_read_var <= nargs
        vn = start_read_var-copts;
       if ischar( varargs{vn} )
           charctr=1;
           if strcmp(varargs{vn}, 'do_avg' )
               doavg = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'avg_groups' )
               avg_groups = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'do_norm' )
               donorm = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'norm_params' )
               params = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'roi_indx' )
               roi_indx = varargs{vn+1};
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
                   roi_indx = varargs{vn};
                   start_read_var = start_read_var+1;
               case 2
                   doavg = varargs{vn};
                   start_read_var = start_read_var+1;
               case 3
                   avg_groups = varargs{vn};
                   start_read_var = start_read_var+1;
               case 4
                   donorm = varargs{vn};
                   start_read_var = start_read_var+1;
               case 5
                   params = varargs{vn};
                   start_read_var = start_read_var+1;
               case 6
                   doplot = varargs{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end

if isempty(avg_groups)
    doavg = false;
end

if donorm
    if isempty(params)
        params.baseline_per_trial=false;
        params.norm_method='percentile';
        params.base_percentile=10;
        params.smooth_type='Moving average';
        params.smooth_scale=100;
        params.acquisition_rate = 10;
    else
        if ~isfield( params, 'baseline_per_trial' )
            params.baseline_per_trial=false;
        end
        if ~isfield( params, 'norm_method' )
            params.norm_method='percentile';
        end
        if ~isfield( params, 'base_percentile' )
            params.base_percentile=5;
        end
        if ~isfield( params, 'smooth_type' )
            params.smooth_type='Moving average';
        end
        if ~isfield( params, 'smooth_scale' )
            params.smooth_scale=100;     
        end
    end
end

end