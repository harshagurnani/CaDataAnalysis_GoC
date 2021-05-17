allvids = struct( 'type', {'ecam','ucam', 'wcam','bcam', 'fcam'}, 'filename', {'EyeCam', 'EyeCam_USB', 'WhiskersCam', 'BodyCam', 'FastCam1'} );
    

for vn = 1:exp.nVid
       M = 0; m=10;
       video = exp.Video.Segments{vn};
       idx = arrayfun( @(jj) ~isempty(strfind( video, allvids(jj).type)),1:5);
       vidtype = allvids(idx).type;
%        for sn=1:exp.nSess
%         sess            = exp.session_types{exp.use_sessions(sn)};
%         m = min([m; VidMI.(sess).(video)(:,1)]);
%         M = prctile([M; VidMI.(sess).(video)(:,1)],99.99);
%        end
       for sn=1:exp.nSess
       sess            = exp.session_types{exp.use_sessions(sn)};
%        VidMI.(sess).(video)(:,1) = (VidMI.(sess).(video)(:,1)-m)/M;
%        VidMI.(sess).(video)(VidMI.(sess).(video)(:,1) > 1,1) = 1;
       VidMI.(sess).(video)(:,2) = VidMI.(sess).(video)(:,2) - videoShift.(sess).(vidtype);
       end
end



%% shift back:

for vn = 1:exp.nVid
       M = 0; m=10;
       video = exp.Video.Segments{vn};
%        idx = arrayfun( @(jj) ~isempty(strfind( video, allvids(jj).type)),1:5);
%        vidtype = allvids(idx).type;

       for sn=1:exp.nSess
       sess            = exp.session_types{exp.use_sessions(sn)};
       VidMI.(sess).(video)(:,2) = VidMI.(sess).(video)(:,2) + videoShift.(sess).(vidtype);
       end
end