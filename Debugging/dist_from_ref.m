function normvector = dist_from_ref( video, ref, TS, varargin)
%% Euclidean Distance from a fixed frame.
% ARGUMENTS :
%   video - Tiff file or array.
%   TS -    Total number of timeframes
%   ref -   Reference frame
% 
% Optional: - Are there split substacks (0 = ImjMacroo, otherise interpreted as number of substacks
%
% RETURNS
%   normvector -    Array of Euclidean distance from reference
%
% Harsha Gurnani. April 2017

normvector = nan(TS,1);

tic
%% If given vector
if isnumeric(video)
        
        s=size(video);
        sr = size(ref);
        FrameDiff = video-ref; 
        if ~( size(s,2)-size(sr,2) == 1 )
            error( 'Reference is not a single time frame OR timeseries array does not have right size. Please check input')
        end
        if size(s,2) == 3
            normvector = squeeze(sum(sum((FrameDiff.*FrameDiff),2),1));
        elseif size(s,2) == 2
            normvector(2:end,1) = squeeze(sum((FrameDiff.*FrameDiff),1));
        elseif size(s,2) == 1
           normvector(2:end,1) = squeeze(sum((FrameDiff.*FrameDiff),1));
        end

%If given file name
elseif ischar(video)
%% Check for substacks
    if nargin == 4
        imj_macro=0;
        n_stacks = varargin{1};

        %%% Parse path and filename
        last_bs = sort(cat(2,find(video== '/'), find(video=='\')));
        if ~isempty(last_bs)
            fpath = video(1:last_bs(end));
            flnm = video(last_bs(end)+1:end);
            if strcmp(flnm(end-3:end), '.tif')
                fltype = 'tif';
                flnm=flnm(1:end-4);
            end
            if strcmp(flnm(end-3:end), '.gif')
                fltype = 'gif';
                flnm=flnm(1:end-4);
            end
        else
            fpath = '';
            flnm = video(1:end);
        end

        if n_stacks == 0
            %%% Used IMJ macro
            imj_macro = 1;
            fls = dir([fpath,flnm, sprintf('-stack*.tiff') ]);
            n_stacks = size(fls,1);
    %         new_flnm = cell(1, n_stacks);
            for jj=1:n_stacks
            new_flnm{jj} = [fpath, flnm, sprintf('-stack%d.tiff', jj)];
            end 
            max_ts = 600;
        else
            new_flnm = cell(1, n_stacks);

            %%% Check filenames with '_'
            fls = dir([fpath, sprintf('*_*%s*.%s',flnm,fltype)]);
            if isempty(fls)
                fls = dir([fpath, sprintf('*%s*_*.%s',flnm,fltype)]);
            end

            for jj=1:n_stacks
                new_flnm{jj} = [fpath, fls(jj).name];
            end
            max_ts = TS/n_stacks;
        end
        

    else 
        n_stacks = 1;
    end

%% Compute distance from reference frame
    if n_stacks > 1

        normvector2 = nan(max_ts, n_stacks);
        parfor jj = 1:n_stacks-1
           normvector2(:,jj) = get_part_dist( new_flnm{jj}, 1, max_ts, ref); 
        end
        normvector = normvector2;
        clear normvector2
        if imj_macro
            last_id = mod(TS, max_ts);
            if last_id == 0
                last_id = max_ts;
            end
            normvector(1:last_id,n_stacks) = get_part_dist( new_flnm{n_stacks}, 1, last_id, ref);
        else
            normvector(:,n_stacks) = get_part_dist( new_flnm{n_stacks}, 1, max_ts, ref); 
        end

        normvector = reshape(normvector, max_ts*n_stacks,1);
        normvector = normvector(1:TS,1);
        
    else
        fltype = video(end-2:end);
        switch fltype
            case 'tif'
                parfor i = 1:TS

                    X1 = imread(video, i);
                    if size(X1,3)==3
                      X1 = rgb2gray(X1);
                    end
                    FrameDiff = squeeze(X1(:,:,1)-ref);
                    normvector(i) = (sum(sum(FrameDiff.*FrameDiff)));%^2;    
                end

        end
    end

    %% Normalize
     M = max(normvector);
     normvector = sqrt(normvector./M); 
     
    toc

end

end

function part_dist = get_part_dist( filenm, start, last, ref)

   part_dist = nan(last-start+1,1);
   for i = start+1:last
        curr_frame = imread(filenm, i);
        if size(ref,3)==3
            curr_frame = rgb2gray(curr_frame);
        end

        FrameDiff = squeeze(curr_frame(:,:,1)-ref);   
        part_dist(i) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

    end

end
