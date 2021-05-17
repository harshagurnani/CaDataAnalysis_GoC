function [ VidMI, Vidtime, ROI ] = set_videoMI_multsess(exp)
% Compute simple motion index for all video segments in all experimental sessions 
    video_split=exp.Video.is_split;
    allvids = struct( 'type', {'ecam','ucam', 'wcam','bcam', 'fcam'}, 'filename', {'EyeCam', 'EyeCam_USB', 'WhiskersCam', 'BodyCam', 'FastCam1'} );
    
    %% which timefile to use
    if ~isfield( exp,'nVid' ),    exp.nVid = numel(exp.Video.Segments);    end
%     if isfield(exp.Video, 'TimeFile')
%         for jj = 1:length(allvids)
%             vidtype = allvids(jj).type;
%             use_vid.(vidtype) =  arrayfun( @(vid) ( (~isempty(strfind( exp.Video.Segments{vid}, vidtype))) || strcmpi( exp.Video.TimeFile{vid}, vidtype) ) , 1:exp.nVid );
%         end
%     else
%         for jj = 1:length(allvids)
%             vidtype = allvids(jj).type;
%             use_vid.(vidtype) =  arrayfun( @(vid)  ~isempty(strfind( exp.Video.Segments{vid}, vidtype)) , 1:exp.nVid );
%         end
%     end
    
    for vid=1:exp.nVid
       for jj = 1:length(allvids)
            vidtype = allvids(jj).type;
            if isfield(exp.Video, 'TimeFile')
                if  (~isempty(strfind( exp.Video.Segments{vid}, vidtype))) || strcmpi( exp.Video.TimeFile{vid}, vidtype), vid_file{vid}=jj; use_vid{vid} = vidtype; end
            else
                if ~isempty(strfind( exp.Video.Segments{vid}, vidtype)), vid_file{vid}=jj; use_vid{vid} = vidtype; end
            end
       end
       if isempty( use_vid{vid} )
           use_vid{vid} = exp.Video.TimeFile{vid}; 
           vid_file{vid} = 0;
           timefile = textread( [exp.Video.(sess), '/', exp.Video.TimeFile{vid}] );
           timefile(:,end) = timefile(:,end) - timefile(1,end);
           Vidtime.(sess).(use_vid{vid}) = timefile;
       end
           
    end
    
    if ~isfield(exp.Video, 'ROI'), exp.Video.ROI = arrayfun( @(jj) 0 , 1:exp.nVid, 'UniformOutput', false); end
    
    %% Process videos
    for sn=1:exp.nSess
        sess            = exp.session_types{exp.use_sessions(sn)}

        % Read all camera timefiles
        for jj = 1:length(allvids)
            vidtype = allvids(jj).type;
            try
            fid = fopen([exp.Video.(sess),'/', allvids(jj).filename,'-relative times.txt' ]);
            tmp = textscan(fid, '%d%d%s', 'delimiter', ' ');
            fclose('all');
            Vidtime.(sess).(vidtype) = double( [tmp{1}, tmp{2}] );
            Vidtime.(sess).(vidtype)(:,2) = Vidtime.(sess).(vidtype)(:,2) - Vidtime.(sess).(vidtype)(1,2); 
            catch
%                 vidtype
            end
        end
         
        
        % compute motion index
        if video_split, folder =  [exp.Video.(sess),'/split/']; 
            % From split substacks saved as *.tif
            for vn = 1:exp.nVid
                video = exp.Video.Segments{vn};
                VidMI.(sess).(video) = simple_motion_index( [folder,video,'.tif'], Vidtime.(sess).(use_vid{vn}), 0);
            end
        else 
            % directly from .avi. Read one frame at a time
            folder = [exp.Video.(sess),'/']; 
            for vn = 1:exp.nVid
                video = exp.Video.Segments{vn};
                MI = Vidtime.(sess).(use_vid{vn});
                if vid_file{vn}~=0, vidFiles = dir( [folder,'*',allvids(vid_file{vn}).filename,'*.avi']);
                else,               vidFiles = dir( [folder,'*',video,'*.avi']);
                end
                nFrames = length(MI);
                c=0;
                for jj = 1:length(vidFiles)                
                   exp.Video.ROI{vn} = get_ROI( [folder,vidFiles(jj).name], exp.Video.ROI{vn} );
%                    [tmpmi, f1, fn, exp.Video.ROI{vn} ] = get_MI_from_video( [folder,vidFiles(jj).name], ROI, nFrames );
%                    oldf1 = f1; oldfn = fn;
%                    MI(c+1:c+length(tmpmi),1) = tmpmi;
%                    if jj>1
%                        d1 = oldfn(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3), :);
%                        d2 = f1(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3), :);
%                        FrameDiff = squeeze(d2(:,:,1)-d1(:,:,1));
%                        MI(c+1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;
%                    end

                    if strcmpi( use_vid{vn}, 'bcam')

                    [tmpmi, ~] = get_MI_from_video( [folder,vidFiles(jj).name], [], false, exp.Video.ROI{vn}, false, false);
                    else
                        % too big (?)
                    [tmpmi, ~] = get_MI_from_video( [folder,vidFiles(jj).name], [], false, exp.Video.ROI{vn}, false, true);
                    end
                    tmpmi=tmpmi{1}(:,1);
                   tmpmi(1) = NaN;
                   MI(c+1:c+length(tmpmi),1) = tmpmi;

                   c = c+length(tmpmi);
                end
                % normalise MI
                 m= 0; %min(MI(MI(:,1)>0,1))
%                 M = max(MI(:,1));
%                 MI(:,1) = (MI(:,1)-m)/M;
                MI = MI(2:end,:);
                VidMI.(sess).(video) = MI;
            end
        end
       
    end
    
    % normalisation across sessions
    for vn = 1:exp.nVid
       M = 0;
       video = exp.Video.Segments{vn};
       for sn=1:exp.nSess
       sess            = exp.session_types{exp.use_sessions(sn)};
%        m = min(VidMI.(sess).(video)(:,1)); 
%        VidMI.(sess).(video)(:,1) = VidMI.(sess).(video)(:,1)-m;
       m = min([m; VidMI.(sess).(video)(:,1)]);
        M = prctile([M; VidMI.(sess).(video)(:,1)],99.99);
       end
       for sn=1:exp.nSess
       sess            = exp.session_types{exp.use_sessions(sn)};
       VidMI.(sess).(video)(:,1) = (VidMI.(sess).(video)(:,1)-m)/M;
       if isfield( exp.Video, 'Timeshift')
           VidMI.(sess).(video)(:,2) = (VidMI.(sess).(video)(:,2)-exp.Video.Timeshift.(sess));
       end
       end
       
    end
    ROI = exp.Video.ROI;
%     
%     disp('.....................Saved Motion Index from Videos')
%   save( [exp.analysed, '\vidMI_', sess, '.mat'], 'VidMI', 'Vidtime','-v7.3' );

end

function ROI = get_ROI( video_file, ROI )

    a = VideoReader( video_file );
    olddata = readFrame(a);
    domore = 'y';
    if isempty(ROI)
        while strcmpi( domore, 'y' ) %%% qqqqq - add more to 
            fig1=figure;imagesc(olddata);
            ROI = floor(getrect(fig1));
            domore = 'n';
        end
        %ROI = [1,1, size(olddata,2)-1, size(olddata,1)-1]; end   
    end
    close all
end


% function [MI, f1, fn, ROI] = get_MI_from_video( video_file, ROI, nFrames )
% 
%     %ROI is a cell array of all ROI from saame video
%     a = VideoReader( video_file );
%     MI  = zeros(nFrames,1);
%     olddata = readFrame(a);
%     f1 = olddata;
%     if isempty(ROI), ROI = [1,1, size(olddata,2)-1, size(olddata,1)-1]; end
%     jj=1;
%     if length(ROI)>1
%         % crop
%         olddata=olddata(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3), :);
%         while hasFrame(a)
%             jj=jj+1;
%             data2 = readFrame(a);
%             data=data2(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3), :);
%             FrameDiff = squeeze(data(:,:,1)-olddata(:,:,1));
%             MI(jj,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;
%             olddata=data;
%         end
%         fn = data2;
%     else
%         % whole frame
%         while hasFrame(a)
%             jj=jj+1;
%             data = readFrame(a);
%             FrameDiff = squeeze(data(:,:,1)-olddata(:,:,1));
%             MI(jj,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;
%             olddata=data;
%         end
%         fn = data;
%     end
%     MI = MI(1:jj);
% 
% end