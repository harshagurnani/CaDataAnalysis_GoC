function [ data, params] = get_rawdata( params )
%% Load raw data either stored in LabView format for Point scans or exported as .mat from .tdms using Antoine's GUI. 
% More flexibility needed for reading other kinds of data.
% NOTE (HG):PARTIAL SUPPORT FOR DATA SAVED WITH ANTOINE's CONTROLLER -
% NEED TO ADD ABILITY TO READ DIRECTLY FROM TDMS FILES - MAYBE USE THE SAME UTILITIES??
% 
% Arguments:
%   - params        Structure with parameters such as 'exp_path'/'data_path', channel, ROIs, trials etc for details about loading
% 
% Returns:
%   - rawdata       3-D for point scans or patch-averaged data. 5-D array for complete patch data. 
%                   3-D format: Dim 1 is Time, Dim 2 is Trial, Dim 3 is ROI(point/patch)
%                   5-D format: Dim 1 is Pixel in Line, Dim 2 is Line Number, Dim 3 is Time, Dim 4 is Trial, Dim 5 is Patchnumber
%   - params        Update params structure with acquisition rate, number of timepoints etc
%

    switch params.import_folder
        
    case 'Labview'
    %% LABVIEW .dat files
        switch params.ROI_type
            case 'poi'
                %% POI
                %Green or red
                if strcmp(params.signal_channel, 'green')
                    channel_let = 'B';
                elseif strcmp(params.signal_channel, 'red')
                    channel_let = 'A';
                end

                %No. of timepoints per trial
                temp = importdata([params.data_path,'ROI-001_POI_Dwell-4.0us_Ch',channel_let,'_Trial-01.dat']);
                params.nTimepoints = size(temp,1);
                clear temp

                data = nan(params.nTimepoints, params.nTrials, params.nROIs );
                %1 page per POI, Columns on each POI page are different trials

                %Channel
                if strcmp(params.signal_channel, 'green')
                    channel_let = 'B';
                elseif strcmp(params.signal_channel, 'red')
                    channel_let = 'A';
                end
                for np = 1:params.nROIs
                    for nt = 1:params.nTrials
                        poi_n = params.ROIs(np);
                        trial_n = params.trials(nt);
                        %%%%%%%%%% need a better way of reading files to
                        %%%%%%%%%% protect against future changes in files
                        %%%%%%%%%% nomenclature 
                        data(:,nt,np) = reshape(importdata( ...
                            [params.data_path,'ROI-',sprintf('%03d',poi_n),'_POI_Dwell-4.0us_Ch',...
                            channel_let,'_Trial-',sprintf('%02d',trial_n),'.dat']), ...
                            [params.nTimepoints,1,1]);
                    end
                end
            case 'patch'
                %%%%%%%%%% READ FROM TDMS _ CURRENTLY NOT SUPPORTED
        end


    otherwise
    %% Used MATLAB exporter: TDMS --> MAT
        %preparing..
%         cd()
        allfiles = dir([params.data_path,'*Scan*']);
        temp = load([params.data_path, allfiles(1).name]);
        
        scan_let = allfiles(1).name(1);
        if strcmp( scan_let,'R')
            scan_type = 'Ribbon';
        elseif strcmp( scan_let, 'S')
            scan_type='Skeleton';
        end
        
        % Time suffix in files can be of variable length for different
        % experiments - read from file directly in case of further
        % format changes
        allUS = find(allfiles(1).name == '_' );
        allUS = allUS(end);
        lastdot = find(allfiles(1).name == '.' );
        TmString = allfiles(1).name(allUS+1:lastdot-1);

        %Channel
        if strcmp(params.signal_channel, 'green')
             channel_let = 2;
        elseif strcmp(params.signal_channel, 'red')
             channel_let = 1;
        end

        switch params.ROI_type
            case 'poi'
                %% POI
               
                %%%%% Have to check the data structure for Skeleton Scan
                params.nTimepoints = size(temp.volume,1);
                clear temp

                data = nan(params.nTimepoints, params.nTrials, params.nROIs );
                %1 page per POI, Columns on each POI page are different trials

                
                for np = 1:params.nROIs
                    for nt = 1:params.nTrials
                        poi_n = params.ROIs(np);
                        trial_n = params.trials(nt);
                        %%%%%%%%%% need a better way of reading files to
                        %%%%%%%%%% protect against future changes in files
                        %%%%%%%%%% nomenclature
                        temp = load( [params.data_path, ...
                                scan_type,'Scan_ROI_', sprintf('%04d',poi_n), '_repeat_' ,...
                                sprintf('%04d',trial_n), '_timepoints_', TmString]); %sprintf('%d',params.nTimepoints)
                        data(:,nt,np) = temp.volume(:,:,:,channel_let);
                    end
                end


            case 'patch'
               %% Patch 

               
               % Size of data matrix for initialisation
               params.nPatchLines = size(temp.volume,2);
               params.nLinePixels = size(temp.volume,1);
               params.nTimepoints = size(temp.volume,3);
               clear temp


               if params.average_patch
                data = nan(params.nTimepoints, params.nTrials, params.nROIs );
               else
                data = nan(params.nLinePixels, params.nPatchLines, params.nTimepoints, params.nTrials, params.nROIs );
                %Each page is one time frame of a patch. Dim 3 = Time series,
                %Dim 4 = Trial, Dim 5 = POI#
               end

               %---------- Note: Original TDMS exported as 1 file per patch per trial
               % More flexibility in future about kind of exports
               for np = 1:params.nROIs
                    for nt = 1:params.nTrials
                        poi_n = params.ROIs(np);
                        trial_n = params.trials(nt);
                        temp = load( [params.data_path, ...
                                scan_type,'Scan_ROI_', sprintf('%04d',poi_n), '_repeat_' ,...
                                sprintf('%04d',trial_n), '_timepoints_', TmString] ); %sprintf('%d',params.nTimepoints)
                        if params.average_patch
                            data(:,nt,np) = squeeze( mean(mean(temp.volume(:,:,:,channel_let))) );
                        else
                            data(:,:,:,nt,np) = temp.volume(:,:,:,channel_let);
                        end
                    end
               end
        end

    end
end