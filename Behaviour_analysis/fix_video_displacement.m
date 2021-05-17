use_ecam = arrayfun( @(vid) ~isempty(strfind( exp.Video.Segments{vid}, 'ecam')), 1:nVid );


for sn=1:nSess

    sess            = session_types{use_sessions(sn)};
    
    % Read times
      
    for vn = 1:nVid
        video = exp.Video.Segments{vn};
        VidMI.(sess).(video)(:,2) = VidMI.(sess).(video)(:,2) - viddisp.(sess).time;
        if use_ecam(vn)
            VidMI.(sess).(video) = VidMI.(sess).(video)(viddisp.(sess).ecam:end,:);
        else
            VidMI.(sess).(video) = VidMI.(sess).(video)(viddisp.(sess).wcam:end,:);
        end
    end

end

save( [exp.analysed, '\vidMI_', sess, '.mat'], 'VidMI', '-v7.3' );
disp('.....................Saved Motion Index from Videos')
