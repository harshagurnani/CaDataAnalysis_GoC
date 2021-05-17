function [xyshift, rho, bestimage ] = MC_with_ECC( RawData, refs, varargin) 

%% Parsing optional arguments

[params, montage] = parse_all_options(varargin, nargin);

%-------------------------------------------------------------------------%
%% Initialising output arrays

nROI = size(RawData, 5);
Tm = size(RawData,3);
nTrials = size(RawData, 4);

if montage
    if strcmp('transform', 'translation')
        xyshift = nan(2, Tm, nTrials);
    end
    bestimage = nan(Tm, nTrials); 
    rho = zeros(Tm, nTrials);    
else
    if strcmp('transform', 'translation')
        all_xyshift = nan(2, Tm, nROI, nTrials);
    end
    bestimage = []; 
    all_rho = zeros(Tm, nROI, nTrials);
end

%-------------------------------------------------------------------------%
%% Calculate optimal transformation
parfor trial = 1:nTrials
    for t=1:Tm
        for roi=1:nROIs
            im = RawData( :,:, t, trial, roi);
            template = refs(:,:,roi);
            [warp, temprho] = iat_ecc(im, template, params);

            if montage
               if temprho > rho(t, trial)
                   rho(t,trial) = temprho;
                   if strcmp('transform', 'translation')
                    xyshift(:, t, trial) = [warp(1,3); warp(2,3)];
                   end
                   bestimage(t,trial) = roi;
               end
            else
                all_rho(t,roi,trial) = temprho;
                if strcmp('transform', 'translation')
                    all_xyshift(:, t, roi, trial) = [warp(1,3); warp(2,3)];
                end
            end
        end
    end
end

if ~montage
    if strcmp('transform', 'translation')
        xyshift = all_xyshift;
    end
    rho = all_rho;
end

end

function [par, montage] = parse_all_options(vargs, nargs)

iterations=[];
levels=[];
transform=[];
initwarp=[];
montage = [];
doMean = [];

if nargs > 2  

    charctr=0;
    numctr=0;
    start_read_var = 3;
    while start_read_var <= nargs
        vn = start_read_var-2;
       if ischar( vargs{vn} )
           charctr=1;
           if strcmp(vargs{vn}, 'iters' )
               iterations = vargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(vargs{vn}, 'montage' )
               montage = vargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(vargs{vn}, 'levels' )
               levels = vargs{vn+1};
               start_read_var = start_read_var + 2; 
           elseif strcmp(vargs{vn}, 'transform' )
               transform = vargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(vargs{vn}, 'initwarp' )
               initwarp = vargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(vargs{vn}, 'do_mean' )
               doMean = vargs{vn+1};
               start_read_var = start_read_var + 2;     
           else 
                error('Argument name %s is not an optional argument', vargs{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   iterations = vargs{vn};
                   start_read_var = start_read_var+1;
               case 2
                   levels = vargs{vn};
                   start_read_var = start_read_var+1;
               case 3
                   initwarp = vargs{vn};
                   start_read_var = start_read_var+1;
               case 4
                   montage = vargs{vn};
                   start_read_var = start_read_var+1;
               case 5
                   doMean = vargs{vn};
                   start_read_var = start_read_var+1;
               case 6
                   transform = vargs{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end


if ~isempty(iterations)
    par.iterations = iterations;
else
    par.iterations = 20;
end

if ~isempty(levels)
    par.levels = levels;
else
    par.levels = 1;
end

if ~isempty(transform)
    par.transform = transform;
else
    par.transform = 'translation';
end

if ~isempty(doMean)
    par.doMean = doMean;
else
    par.doMean=false;
end

if ~isempty(initwarp)
    par.initwarp = initwarp;
end

if isempty(montage)
    montage=true;
end

end