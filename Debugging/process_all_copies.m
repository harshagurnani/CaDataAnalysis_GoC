function [last_ts, dups, MI, refdist, crosscorr, corrlag, speed_corr, speed_lag ] = process_all_copies( indx, vd_path, tiffiles, timefiles, partlabel, varargin )
%%    Arguments:
%
%       - indx :        Row vector with values from {1,2,3,...,N}
%                       corresponding to reading the right tif and time
%                       files. Currently N <= 3
%       - vd_path:      Path for all videos (tif AND time files)
%       - tiffiles:     Cell array with tif file names/path
%       - timefiles:    Cell array with time file names/path
%       - partlabel:    Label of body part - from 'wheel', whiskpad', extra1', 'extra2'
%       
%   Optional:
%       - speed data:   Only give for wheel otherwise correlations will be
%       calculated between tiffile data and speed
%
%   RETURNS:
%   
%       - last_ts:      Cell array with LAST TIMESTAMP for each video -
%                       indexed as {video_indx}
%       - dups:         Structure to store duplications in timestamp or
%                       video frame
%                  dups.timestamps - Nx2 Array with NUMBER OF duplications
%                       in each timestamp file in columns 1(frame no.) and
%                       2 (actual timestamp)  -order could be reversed in
%                       older files. Indexed as: (video_indx, column_num)
%                  dups.frames - Cell array - Frame number/index were the
%                       motion index for the corresponding video is 0
%                       (normed MI = 0 suggests possible duplicated frames
%                       - especially if there is a block of such frames,
%                       particularly in later recordings). 
%                       Indexed as: {video_indx}
%       - MI:           Cell array with (T-1)x2 "MI array" for each video -
%                       Column 1 of each video's MI array is motion index,
%                       column 2 is frame timestamp. 
%                       Indexed as: {video_indx}( frame_num, 1 or 2)
%       - RefDist:      Cell array with (T-1)x2 "Distance from Ref array"
%                       for each video. Column 1 of each video's RefDist
%                       array is Euclidean dist from first frame
%                       (reference), column 2 is frame timestamp. Indexed
%                       as: {video_indx}( frame_num, 1 or 2)
%       - crosscorr:    Cell-array for lag cross-correlations between
%                       different videos (if N>1) - indexed as
%                       {video_indx1}{video_indx2}
%       - corrlag:      Cell-array for lag values corresponding to crosscorr 
%                       between different videos (if N>1) - indexed as
%                       {video_indx1}{video_indx2}
%       - speed_corr:   Cell array - Each cell has lag correlations between
%                       speed and correspondingly indexed video
%                       Indexed as: {video_indx}
%       - speed_lag:    Cell array - Each cell has lag values at which
%                       correlations were calculated between speed and
%                       correspondingly indexed video
%                       Indexed as: {video_indx}
%-------------------------------------------------------------------------%
    nFiles = size(indx,2);
    if nargin > 5
        speed = varargin{1};
    else
        speed = [];
        
    end
    speed_corr={};
    speed_lag = {};
    
    last_ts={}; dups={}; MI={}; refdist = {}; crosscorr={}; corrlag={};
    
    if nFiles>=1

        MI = cell(1, nFiles);
        refdist = cell(1, nFiles);
        plot_mi = {};

        for jj=1:nFiles
           
           %if given as -  path/file;0
           if ~isempty(find(tiffiles{jj}==';', 1))
               semic = find(tiffiles{jj}==';', 1);
               
               nstck = str2double(tiffiles{jj}(semic+1:end));
               tiffiles{jj} = tiffiles{jj}(1:semic-1);
               
           else 
               nstck = [];
           end
           
           time_semic = find(timefiles{jj}==';',1);
           if ~isempty(time_semic)
                   timefiles{jj} = timefiles{jj}(1:time_semic-1);
           end
           
           % Parse standard format for tif data
           if size(tiffiles{jj},2) <= 4
               if size(tiffiles{jj},2) > 1
                nstck = str2double(tiffiles{jj}(2:end));  
               end
               switch tiffiles{jj}(1)
                   case 'e'
                       tiffiles{jj} = sprintf('EyeCam_%s.tif', partlabel);
                   case 'w'
                       tiffiles{jj} = sprintf('WhiskersCam_%s.tif', partlabel);
                   case 'b'
                      tiffiles{jj} = sprintf('BodyCam_%s.tif', partlabel);
               end
           end
           
           %Parse standard format for time data
           if size(timefiles{jj},2) <= 4
               if size(timefiles{jj},2) > 1 && isempty(nstck)
                nstck = str2double(timefiles{jj}(2:end));  
               end
               switch timefiles{jj}(1)
                   case 'e'
                       timefiles{jj} = 'EyeCam-relative times.txt';
                   case 'w'
                       timefiles{jj} = 'WhiskersCam-relative times.txt';
                   case 'b'
                       timefiles{jj} = 'BodyCam-relative times.txt';
               end
           end

           datatimes{jj} = importdata([vd_path, timefiles{jj}]);
           last_ts{jj} = datatimes{jj}(end,2);

           % Number of Duplications in timestamps
           dups.timestamps(jj, :) = [size(datatimes{jj}(:,1),1)-size(unique(datatimes{jj}(:,1)),1), ... no. of duplications in frame number
                                    size(datatimes{jj}(:,2),1)-size(unique(datatimes{jj}(:,2)),1)]; %   no. of duplications in timestamp

           disp(sprintf('Computing motion index for %s %d ...', partlabel, jj)) %#ok<*DSPS>
           
           try 
               if ~isempty(nstck)
                        MI{jj} = simple_motion_index( [vd_path, tiffiles{jj}], datatimes{jj} , nstck);
               else 
                   MI{jj} = simple_motion_index( [vd_path,tiffiles{jj}], datatimes{jj} );
               end

               % Possible duplicated frames: Indices where MI for wheel video jj is zero. 
               dups.frames{jj} = find(MI{jj}(:,1) == 0 );   
           catch ME
                disp(sprintf('%s %d MI: %s',partlabel, jj,ME.message))
           end
           disp(sprintf('Computing distance from reference frame for %s %d ...', partlabel, jj)) 
           tic
           %% UNCOMMENT IF YOU NEED DISTANCE FROM FIRST FRAME - eg for periodic movements
%            refframe = imread( [vd_path, tiffiles{jj}], 1);
%            refframe = squeeze(refframe);
%            refdist{jj} = dist_from_ref( [vd_path, tiffiles{jj}], refframe, datatimes{jj}(end,1) );

           toc
           clear nstck
        end

        if ~isempty(speed)
            disp(sprintf('\n Computing correlations between %s motion as well as speed ...', partlabel ))
        else
            disp('No speed data available. Cannot compute correlations against speed.');
        end
    else
        disp(sprintf('\n No %s data available.',partlabel ))
    end


    %% Correlations for wheel and speed data
    % ----------------------------------------- 

    ds_rate = 25;   %Hz - can be higher, if cameras recorded at higher frame rate.
    lagbins = 1000/ds_rate;   %No. of lag bins on each side of 0


    plot_ws={};     % Indices of successful speed-tiffile corr
    plot_ww={};     % Indices of tiffiles for cross corr

    if nFiles>0 
    for jj=1:size(tiffiles,2)

        if  ~isempty(speed)
            try
            [speed_corr{jj}, speed_lag{jj}] = get_lagcorr_ds_synctest( MI{jj}(:,1), MI{jj}(:,2),...
                                                                abs(speed(:,2)), speed(:,1), ds_rate, lagbins ) ;                                                 
            plot_ws = cat(2, plot_ws, num2str(jj));
            catch ME
                disp(sprintf('Error with correlation between speed and %s %d motion index:\n %s', partlabel, jj, ME.message))
            end
        end

        

        for kk=jj+1:size(tiffiles,2)
            try
           [ crosscorr{jj}{kk}, corrlag{jj}] = get_lagcorr_ds_synctest( MI{jj}(:,1), MI{jj}(:,2),...
                                                                    MI{kk}(:,1), MI{kk}(:,2),...
                                                                    ds_rate, lagbins); 

            plot_ww = cat(2, plot_ww, sprintf('%s %d - %s %d', partlabel, jj, partlabel, kk) );
            catch ME
                disp(sprintf('Error with correlation between %s %d and %s %d motion index:\n %s', partlabel, jj, partlabel, kk, ME.message))
            end
        end
    end
    end



    %% Plot correlation for speed and wheel data
    % ----------------------------------------- 

    % 1) Lag correlation between speed data and motion index of wheel
    if nFiles>=1
        if ~isempty(speed)
            figure;colormap(lines)
            title(sprintf('Lag correlation between speed data and motion index of %s',partlabel))
            for jj=1:nFiles
                try
                plot(speed_lag{jj}, speed_corr{jj})
                catch ME
                end
                hold on
            end
            legend(strcat(partlabel, plot_ws))
            xlabel( 'Lag (ms)')
            ylabel('Correlation coefficient')
        end
        
        figure;colormap(lines)
        title(sprintf('Speed data and motion index of %s', partlabel))
        if ~isempty(speed)
            plot( speed(:,1), 4*smooth(abs(speed(:,2)),30)/max(abs(speed(:,2)))+0.1 )
            lg = cat(2, 'Speed' ,strcat(partlabel,plot_mi));
        else
            lg = strcat(partlabel,plot_mi);
        end
        hold on
        for jj=1:nFiles
            try
            plot(MI{jj}(:,2), MI{jj}(:,1))
            catch ME
            end
            hold on
        end
        legend(lg)
        xlabel( 'Time (ms)')
        ylabel('Normalised speed or motion index')
    end

    % 2) Lag correlation between motion indices of different wheel videos
    if nFiles>1
        figure;colormap(lines)
        title(sprintf('Lag correlation between motion indices of different %s videos',partlabel))
        for jj=1:nFiles
            for kk=jj+1:nFiles
                try
                plot(corrlag{jj}, crosscorr{jj}{kk})
                catch ME
                end
                hold on
            end
        end
        legend(plot_ww)
        xlabel( 'Lag (ms)')
        ylabel('Correlation coefficient')
    end

end