function Speed = run_speed_multsess( exp )
    for sn = 1:exp.nSess
        sess            = exp.session_types{exp.use_sessions(sn)};
        Speed.(sess)    = cell2mat( (get_speed_data([exp.(sess), '\Speed_Data\']))' );
    end
    save( [exp.analysed, '\Speed.mat'], 'Speed', '-v7.3' );
    disp('.....................Saved Speed Data')
end