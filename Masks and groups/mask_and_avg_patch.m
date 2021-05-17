function [ masked_data, time_mask, masks ] = mask_and_avg_patch( data, times, cropping, varargin)
%% Create and appy masks to patches, and get average traces of masked patch. Normalisation is optional.
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
nROI = size(data,5);

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

%% Parse optional arguments
[ fthresh, dil_factor, donorm, params, doplot ] = parse_options( nargin, varargin, 3 );

%% Create masks based on average patch image
im={};
for roi=1:nROI
    im{roi} = mean( reshape(data(:,:,:,:,roi), [ny,nx, Tm*nTrials]), 3 );
end

masks={};
for jj=1:nROI
    [~,threshold]=edge(im{jj},'sobel');     
    BW=edge(im{jj}, 'sobel',fthresh*threshold);         % Detect edges
    BWdil = imdilate(BW, strel('square',dil_factor));   % Smoothen/expand edges
    masks{jj} = imfill(BWdil, 'holes');                 % fill in any outlined objects
end


%% Apply masks and average
for roi=1:nROI
    temp = reshape(data(ycrop,xcrop,:,:,roi), [ny*nx, Tm, nTrials]);    
    allx = floor(find( masks{roi} == 1 )/ny)+1;
    ally = mod(find( masks{roi} == 1 ),ny)+1;
    
    c_time = [floor(mean(ally)) floor(mean(allx))];
    

    for trial=1:nTrials
    for t=1:Tm
        masked_data(t,trial,roi) = mean( temp(masks{jj}(:), t, trial), 1 );
        time_mask(t,trial,roi) = times(ycrop(c_time(1)), xcrop(c_time(2)), t, trial, roi);
    end
    end
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

function [ fthresh, dil_factor, donorm, params, doplot ] = parse_options( nargs, varargs, copts)

fthresh=[];
dil_factor = [];
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
           if strcmp(varargs{vn}, 'fthresh' )
               fthresh = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'dilate' )
               dil_factor = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'do_norm' )
               donorm = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'norm_params' )
               params = varargs{vn+1};
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
                   donorm = varargs{vn};
                   start_read_var = start_read_var+1;
               case 2
                   params = varargs{vn};
                   start_read_var = start_read_var+1;
               case 3
                   doplot = varargs{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end

if isempty(fthresh)
    fthresh = 0.6;
end

if isempty(dil_factor)
    dil_factor=3;
end

if isempty(donorm)
    donorm=false;
end

if isempty(doplot)
    doplot=false;
end

if donorm
    if isempty(params)
        params.baseline_per_trial=false;
        params.norm_method='percentile';
        params.base_percentile=5;
        params.smooth_type='Moving average';
        params.smooth_scale=100;
    else
        if ~isfield( params, 'baseline_per_trial' )
            params.baseline_per_trial=false;
        end
        if ~isfield( params, 'norm_method' )
            params.norm_method='percentile';
        end
        if ~isfield( params, 'base_percentile' )
            params.base_percentile=10;
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