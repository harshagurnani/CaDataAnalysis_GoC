for sn = 1:nSess
    sess = session_types{use_sessions(sn)};
    a = ImagingData( exp.(sess),[]);
    a.add_itt;
    a1 = ImagingData_P1( exp.(sess),a, []);
    a1.add_itt;
    a1.add_encoder_data;
    Speed.(sess) = a1.speed;
end
    
    
    