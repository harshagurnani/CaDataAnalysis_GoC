function [ VidMI, Vidtime ] = run_videoMI_multsess(exp)
% Compute simple motion index for all video segments in all experimental sessions 
    video_split=true;

    % which timefile to use
    if ~isfield( exp,'nVid' ),    exp.nVid = numel(exp.Video.Segments);    end
    if isfield(exp.Video, 'TimeFile')
        use_ecam = arrayfun( @(vid) ( (~isempty(strfind( exp.Video.Segments{vid}, 'ecam'))) || strcmpi( exp.Video.TimeFile{vid}, 'ecam') ) , 1:exp.nVid );
        use_wcam = arrayfun( @(vid) ( (~isempty(strfind( exp.Video.Segments{vid}, 'wcam'))) || strcmpi( exp.Video.TimeFile{vid}, 'wcam') ) , 1:exp.nVid );
        use_bcam = arrayfun( @(vid) ( (~isempty(strfind( exp.Video.Segments{vid}, 'bcam'))) || strcmpi( exp.Video.TimeFile{vid}, 'bcam') ) , 1:exp.nVid );
    else
        use_ecam = arrayfun( @(vid)  ~isempty(strfind( exp.Video.Segments{vid}, 'ecam')) , 1:exp.nVid );
        use_wcam = arrayfun( @(vid)  ~isempty(strfind( exp.Video.Segments{vid}, 'wcam')) , 1:exp.nVid );
        use_bcam = arrayfun( @(vid)  ~isempty(strfind( exp.Video.Segments{vid}, 'bcam')) , 1:exp.nVid );

    end
    for sn=1:exp.nSess

        sess            = exp.session_types{exp.use_sessions(sn)};

        % Read all camera timefiles
        try
        Vidtime.(sess).ecam = textread( [exp.Video.(sess),'\EyeCam-relative times.txt' ] );
        Vidtime.(sess).ecam(:,2) = Vidtime.(sess).ecam(:,2) - Vidtime.(sess).ecam(1,2); 
        catch
        end
        try
        Vidtime.(sess).wcam = textread( [exp.Video.(sess),'\WhiskersCam-relative times.txt' ]);
        Vidtime.(sess).wcam(:,2) = Vidtime.(sess).wcam(:,2) - Vidtime.(sess).wcam(1,2);
        catch
        end
        try
        Vidtime.(sess).bcam = textread( [exp.Video.(sess),'\BodyCam-relative times.txt' ]);
        Vidtime.(sess).bcam(:,2) = Vidtime.(sess).bcam(:,2) - Vidtime.(sess).bcam(1,2);
        catch
        end
        
        % compute motion index
        for vn = 1:exp.nVid
            video = exp.Video.Segments{vn};

            if video_split, folder =  [exp.Video.(sess),'\split\']; else, folder = [exp.Video.(sess),'\']; end
            if use_ecam(vn)
                VidMI.(sess).(video) = simple_motion_index( [folder,video,'.tif'], Vidtime.(sess).ecam, 0);
            elseif use_wcam(vn)
                VidMI.(sess).(video) = simple_motion_index( [folder,video], Vidtime.(sess).wcam, 0);
            elseif use_bcam(vn)
                VidMI.(sess).(video) = simple_motion_index( [folder,video], Vidtime.(sess).bcam, 0);
            else
                timefile = textread( [exp.Video.(sess), '\', exp.Video.TimeFile{vn}] );
                timefile(:,end) = timefile(:,end) - timefile(1,end);
                VidMI.(sess).(video) = simple_motion_index( [folder,video], timefile, 0);
                Vidtime.(sess).other.(video) = timefile;
                [~,ia,~] = unique(VidMI.(sess).(video)(:,2));
                VidMI.(sess).(video) = VidMI.(sess).(video)(ia,:);

            end
        end

    end

    save( [exp.analysed, '\vidMI_', sess, '.mat'], 'VidMI', 'Vidtime','-v7.3' );
    disp('.....................Saved Motion Index from Videos')
end