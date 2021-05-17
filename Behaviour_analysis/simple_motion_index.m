function [MI] = simple_motion_index(video, Timestamp, varargin)
%% Determine a motion index based on Jelitai et al., Nature communications,
% 2016.
% Note. HG 1/3/2017: Adapted from Fred's code.
%
% the video is a tif file (grayscale). The time stamps of the video are used to
% determine the number of frame and are appended to the matrix MI. Timestamp
% should be a one column matrix. Optional argument can be given if original
% tiff has been split into substacks.
%
% Usage: 
%---------------------
%     --- A     Reading video from a tif file
%     --------------------------------------------------------------------
%       
%         [MI] = simple_motion_index('example.tif', time)      reads a single file specified by the complete path/name, 
%                                                              uses the timestamps in the corresponding time vector
%               OR
%
%     --- B     With file split into substack files (parallelizes reading)
%     --------------------------------------------------------------------           
%           - Give additional third argument which decides how substacks
%           were created, and how many of them to load.
%           
%         
%    -   [MI] = simple_motion_index('example.tif', time, 0 )    if tiff has been split using my ImjJ Macro
%                                                              so the naming of stacks is as I wrote it
%               OR
%
%
%    -   [MI] = simple_motion_index('example.tif', time, -1)   if tiff has been split into substacks with same prefix/suffix as video
%                                                               - read ALL such files
%                                                              Files are identified by looking for the character '_' or '-' in the filename as well as the filename part of the given video. 
%                                                              eg. example_subset1.tif, example_subset2.tif etc 
%                                                         OR   eg. stack1_example.tif,  stack2_example.tif etc
%
%               OR
%    -   [MI] = simple_motion_index('example.tif', time, K )   if tiff has been split into 'K' substack files (or only 'K' files need to be read) with same prefix/suffix as video
%                                                              Same as above.                                                                    
%                                                                        
%               OR
%
%   --- C       Where video is an array in the workspace
%   ---------------------------------------------------------------------
%        [MI] = simple_motion_index(frames, time)              where frames is an array (XxYxT, or XxT, or Tx1 (Values of some variable in time) )
%
%   ---------------------------------------------------------------------
%   Useful Tips
%   ---------------------------------------------------------------------
%       - if there are similar files but you only want to process one of
%       them, give the exact filename WITH EXTENSION as the first argument
%       Eg. Whiskers1.tif 
%           or Whiskers1 ( if no other file has the fragment Whiskers1)
%
%       - if there are multiple whiskers video that need to be
%       CONCATENATED, give the common partial filename and NO EXTENSION. 
%       For example, for combining Whiskers1 and Whiskers2, give the first argument as Whisk
%       or Whisker etc.. (depending on other files in the directory)
%       The files will be sorted, so all stacks of Whiskers1_ ...  will be
%       followed by sorted substacks of Whiskers2_ ... (even if each
%       partitioned video file had relative numbering).
%
%       - However, DO NOT use this file to combine two parallelly acquired
%       videos (Eg: two cropped sections from the same video). Parfor is
%       already running as many threads as you have asked. Calling this
%       function twice on two separate files does not slow you down.
%
%       - The number of timepoints kept will be the minimum of those in the
%       timestamp vector, or the total number of frames in specified number
%       of stacks. If you send a cropped time vector, the mi will then be
%       chopped. However, this chopping is done at the end to accomodate
%       for stacks of variable lengths.
%
%       - The stacks are no longer assumed to have same number of frames.
%       
%
% ARGUMENTS :
%       video    -     Can be (path+)filename or a MATLAB array (Last dim is interpreted as time)
%                      If extension is given, it is assumed to be an exact
%                      filename, otherwise we look for partial filename
%                      match.
%       Timestamp-     (1 or) 2-column array of Frame number and timestamps
%
% Optional:
%       Varargin{1} - If 0, my Macro naming is assumed. If -1, all matching files are processed.
%                     Otherwise interpreted as the number of files with same prefix as video to be
%                     read as successive substacks.
%
% RETURNS:
%   MI -    2 column array. First column is motion index (between 0 and 1).
%           Second column are timestamps copied from given Timestamp vector.
%
%
%
% Harsha Gurnani. March 2017

nFrames = length(Timestamp);



tic
%% If given vector
if isnumeric(video)
        
        s=size(video);
        nFrames_to_keep = s(end)-1;
        if size(s,2) == 3
            FrameDiff = video(:,:,2:end)-video(:,:,1:end-1);
            MI(2:end,1) = squeeze(sum(sum((FrameDiff.*FrameDiff),2),1));%^2;
        elseif size(s,2) == 2
           FrameDiff = video(:,2:end)-video(:,1:end-1);
            MI(2:end,1) = squeeze(sum((FrameDiff.*FrameDiff),1));%^2;
            
        elseif size(s,2) == 1
           FrameDiff = video(2:end)-video(1:end-1);
           MI(2:end,1) = FrameDiff.*FrameDiff;%^2;
        end

%% If given file name
elseif ischar(video)
% Check for substacks
    if nargin == 3
        imj_macro=0;
        n_stacks = varargin{1};

        %%% Parse path and filename
        last_bs = sort(cat(2,find(video== '/'), find(video=='\')));
        if ~isempty(last_bs)
            fpath = video(1:last_bs(end));
            flnm = video(last_bs(end)+1:end);
        else
            fpath = '';
            flnm = video(1:end);
            
        end
        if strcmp(flnm(end-3:end), '.tif')
            fltype = 'tif';
            flnm=flnm(1:end-4);
            partial = false;
        elseif strcmp(flnm(end-4:end), '.tiff')
            fltype = 'tif';
            flnm=flnm(1:end-5);
            partial = false;
        elseif strcmp(flnm(end-3:end), '.gif')
            fltype = 'gif';
            flnm=flnm(1:end-4);
            partial = false;
        end
        if ~exist('fltype', 'var'),     partial = true; fltype = '.tif';    end     %If no filetype mentioned, tif(f) by default, and assumed to be partial filename.
        
        if n_stacks == 0
            %%% Used IMJ macro
            imj_macro = 1;
            if partial
                fls = dir([fpath, '*', flnm, sprintf('*-stack*.tiff') ]);   %to allow for partial match
            else
                fls = dir([fpath, flnm, sprintf('-stack*.tiff') ]);         %exact filename given
            end
            n_stacks = size(fls,1);
            
        else
            %%% Substacks created independently
            %%% Check filenames with '_' - filename (flnm) needs to be
            %%% consistently prefix or suffix.
            if partial
                fls = dir([fpath sprintf('*_*%s*.%s*',flnm,fltype)]);
                if isempty(fls)
                    fls = dir([fpath sprintf('*%s*_*.%s*',flnm,fltype)]);
                end
                if isempty(fls)
                    fls = dir([fpath sprintf('*%s*-*.%s*',flnm,fltype)]);
                end
                if isempty(fls)
                    fls = dir([fpath sprintf('*%s*-*.%s*',flnm,fltype)]);
                end
            else
                %Exact filename - only *NNN_/-FILENAME.EXT or
                %FILENAME_/-*NNN.EXT will be read
                fls = dir([fpath sprintf('*_%s.%s*',flnm,fltype)]);
                if isempty(fls)
                    fls = dir([fpath sprintf('%s_*.%s*',flnm,fltype)]);
                end
                if isempty(fls)
                    fls = dir([fpath sprintf('*-%s.%s*',flnm,fltype)]);
                end
                if isempty(fls)
                    fls = dir([fpath sprintf('%s-*.%s*',flnm,fltype)]);
                end                
            end
            
            if n_stacks == -1,   n_stacks = numel(fls);   end       %Read all files?
            
            
        end
        
        % sort files 
        % (successive video names need to follow dictionary order - so either stack numbers are serial, or they are used as suffix!)
        allfilenames = string( {fls.name} );
        [~,idx] = sort( allfilenames );
        fls = fls(idx);

        new_flnm = cell(1, n_stacks);
        for jj=1:n_stacks
            new_flnm{jj} = [fpath fls(jj).name];
        end

        stack_nFrames = arrayfun( @(jj) numel( imfinfo([fpath fls(jj).name]) ), 1:n_stacks);      %No. of timeframes in each substack.
        max_ts = max( stack_nFrames); total_ts = sum( stack_nFrames);
        nFrames_to_keep = min(nFrames, total_ts );
    else 
        n_stacks = 1;
        nFrames_to_keep = nFrames;
    end

    MI = nan(nFrames_to_keep,2);    
    %% Compute motion index
    if n_stacks > 1

        partial_mi = cell(n_stacks,1);
        parfor jj = 1:n_stacks
           partial_mi{jj} = get_part_MI( new_flnm{jj}, 1, stack_nFrames(jj)); 
           first_frame{jj} = imread(new_flnm{jj},1);
           last_frame{jj} = imread(new_flnm{jj},stack_nFrames(jj));
        end

%         last_id = mod(nFrames, max_ts);
%         if last_id == 0 || nFrames_to_keep<nFrames
%             last_id = max_ts;
%         end
%         
% %         if imj_macro            
%         partial_mi(1:last_id,n_stacks) = get_part_MI( new_flnm{n_stacks}, 1, last_id);
%         first_frame{n_stacks} = imread(new_flnm{n_stacks},1);
%         last_frame{n_stacks} = imread(new_flnm{n_stacks},last_id);
        
%         else
%             partial_mi(:,n_stacks) = get_part_MI( new_flnm{n_stacks}, 1, max_ts); 
%             first_frame{n_stacks} = imread(new_flnm{n_stacks},1);
%             last_frame{n_stacks} = imread(new_flnm{n_stacks},max_ts);
%         end

       
        %%% Compute MI at junctions of stacks
        for jj=2:n_stacks
           FrameDiff = squeeze(first_frame{jj}(:,:,1)-last_frame{jj-1}(:,:,1));
           partial_mi{jj}(1) = (sum(sum(FrameDiff.*FrameDiff)));
        end
        partial_mi = cell2mat(partial_mi);%reshape(partial_mi, max_ts*n_stacks,1);

        MI(:,1) = partial_mi(1:nFrames_to_keep);

    else
        %% Warning - this part has not been tested for bugs with options yet. Use with caution.
        deci = strfind( video, '.');
        if ~isempty( deci),
        fltype = video(deci+1:deci+3);
        else
            flytpe = '.tif';
        end
        
        switch fltype
            case 'tif'
                %%% Normal case with single tif file
                parfor i = 2:nFrames_to_keep

                    X = imread(video,i-1);
                    X1 = imread(video, i);

                    FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));

                    MI(i,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                end

            case 'gif'
                %%% Too slow for loading in memory. DO NOT USE
                nframes = 100;
                start=1;
                ncycles = ceil(size(Timestamp,1)/nframes);
                for cyc=1:ncycles
                    last = min(start+nframes, size(Timestamp,1));
                    curr_stack = start:last;
                    vid_part = imread(video, curr_stack);
                    if cyc==1
                        temp = nan(nframes,1);
                        parfor i=2:nframes

                            FrameDiff = squeeze(vid_part(:,:,1,i-1)-vid_part(:,:,1,i));

                            temp(i,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                        end
                    else  
                        temp = nan(nframes,1);
                        parfor i=2:nframes+1

                            FrameDiff = squeeze(vid_part(:,:,1,i-1)-vid_part(:,:,1,i));

                            temp(i-1,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                        end

                    end
                    MI(curr_stack(1:end-1),1) = temp;
                    start=last;
                end

        end
    end


end

%% Normalize

m = 0;%min(MI(1:end-1000,1));
M = max(MI(:,1));

MI(:,1) = (MI(:,1)-m) ./ (M-m); 

toc

MI(:,2) = Timestamp(1:nFrames_to_keep,end);
MI = MI(2:end,:);

end

function part_mi = get_part_MI( filenm, start, last)

   part_mi = nan(last-start+1,1);
   for i = start+1:last

        X = imread(filenm,i-1);
        X1 = imread(filenm, i);

        FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));   
        part_mi(i) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

    end

end