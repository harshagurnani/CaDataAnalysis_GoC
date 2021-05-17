function [ masked_data,  varargout ] = mask_nonbinary_and_avg( data, times, masks, varargin)
%% Apply nonbmasks to patches given by roi numbers (so you can use multiple masks in same patch). Normalisation and averaging groups is optional.
% 
% ARGUMENTS:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%

nx = size(data,1);
ny = size(data,2);
%reshape into 4D: x-y-t-roi
if length(size(data))==5, data = reshape(data, [nx,ny,size(data,3)*size(data,4),size(data,5)]); end
Tm = size(data,3);
nr = size(data,4);

if ~isempty(times) , do_time  = true; else, do_time = false; end     %Crop times matrix?

%% Parse optional arguments
[ roi_indx, cropping, doavg, avg_groups, donorm, params, doplot, dobinary ] = parse_options( nr, varargin{:} );
nROI = length(roi_indx) ;

masked_data = nan( Tm, nROI);
if do_time, masked_time = nan( Tm, nROI);  end

if size(cropping,2) == 1
    DX = cropping;
    DY = cropping;
else
    DX = cropping(1);
    DY = cropping(2);
end

xcrop = DX+1:nx-DX;
ycrop = DY+1:ny-DY;

nx = size(xcrop,2);
ny = size(ycrop,2);



%% Apply masks
if do_time
    % crop data and timestamps
    for roi=1:nROI
    patch_id = roi_indx(roi);
    %Time only uses pixels with positive weights (negative weights don't count towards signal but background subtraction )
    if dobinary, masks{roi} = (masks{roi}>0);end
    time_mask = masks{roi}.*(masks{roi}>0);   
    nPix = sum(time_mask(:));
    for jj =1:Tm
        masked_data(jj,roi) = nansum(nansum( masks{roi}.*data(xcrop,ycrop,jj,patch_id)));
        masked_time(jj,roi) = nansum(nansum( time_mask.*times(xcrop,ycrop,jj,patch_id)))/nPix;   %Averaged time
    end
    end
else
    % crop only data
    for roi=1:nROI
    patch_id = roi_indx(roi);
    for jj =1:Tm
        masked_data(jj,roi) = nansum(nansum( masks{roi}.*data(xcrop,ycrop,jj,patch_id)));
    end
    end    
end

%% Optional averaging
if doavg
    ngroups = numel(avg_groups);
    avg_data = nan( Tm,  ngroups);
    if do_time
        avg_time = nan( Tm,  ngroups);   
        for group=1:ngroups
            % Mean or sum?
            avg_data(:,group) = mean( masked_data(:, avg_groups{group}),3 );
            avg_time(:,group) = mean( masked_time(:, avg_groups{group}),3 );
        end
        masked_data = avg_data;
        masked_time = avg_time;
    else
        for group=1:ngroups
            % Mean or sum?
            avg_data(:,group) = mean( masked_data(:, avg_groups{group}),3 );
        end
        masked_data = avg_data;
    end
    nROI = ngroups;
end

if nargout>1
    varargout{1} = masked_time;
end

%% Optional normalisation and plotting
if donorm
    if do_time
        [norm, ~] =normalize_and_smooth( reshape(masked_data, [Tm,1,nROI]), reshape(masked_time, [Tm,1,nROI]), params);
        norm = reshape(norm, [Tm, nROI]);
        if doplot
            plot(masked_time, norm)
        end
        masked_data = norm;
    else
        warning('No timestamps given, cannot normalize');
    end

end

end





%-------------------------------------------------------------------------%
%% PARSING OPTIONS %%
%-------------------------------------------------------------------------%

function [ roi_indx, cropping, doavg, avg_groups, donorm, params, doplot, dobinary ] = parse_options( nroi,varargin )

cropping = 0;
roi_indx = 1:nroi;
doavg=false;
avg_groups = [];
donorm = false;
params = [];
doplot =false;
dobinary = false;

nvarargs = length(varargin);
if nvarargs>0  

    charctr=0;
    numctr=0;
    vn = 1;
    while vn <= nvarargs
       if ischar( varargin{vn} )
           charctr=1;
           if strcmp(varargin{vn}, 'do_avg' )
               doavg = varargin{vn+1};
               vn = vn + 2;
           elseif strcmp(varargin{vn}, 'avg_groups' )
               avg_groups = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'do_norm' )
               donorm = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'norm_params' )
               params = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'roi_indx' )
               roi_indx = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'do_plot' )
               doplot = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'cropping' )
               cropping = varargin{vn+1};
               vn = vn + 2;  
           elseif strcmp(varargin{vn}, 'do_binary' )
               dobinary = varargin{vn+1};
               vn = vn + 2;  
               
           else 
                error('Argument name %s is not an optional argument', varargin{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   roi_indx = varargin{vn};
                   vn = vn+1;
               case 2
                   doavg = varargin{vn};
                   vn = vn+1;
               case 3
                   avg_groups = varargin{vn};
                   vn = vn+1;
               case 4
                   donorm = varargin{vn};
                   vn = vn+1;
               case 5
                   params = varargin{vn};
                   vn = vn+1;
               case 6
                   doplot = varargin{vn};
                   vn = vn+1;
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