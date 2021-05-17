beh_all.rwristX = [];
beh_all.rwristY = [];

beh_all.vel_rwristX = [];
beh_all.vel_rwristY = [];


interSess = 1000;
maxt=0;

for jj=1:nSess
sess = session_types{use_sessions(jj)};
totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;

% concatenate behaviour
beh_all.rwristX = [ beh_all.rwristX ; Limb.(sess).RWristX+ [maxt-totalt, 0] ];
beh_all.rwristY = [ beh_all.rwristY ; Limb.(sess).RWristY+ [maxt-totalt, 0] ];

beh_all.vel_rwristX = [ beh_all.vel_rwristX ; [Limb.(sess).RWristX(2:end,1)+ maxt-totalt, diff(Limb.(sess).RWristX(:,2))] ];
beh_all.vel_rwristY = [ beh_all.vel_rwristY ;  [Limb.(sess).RWristY(2:end,1)+ maxt-totalt, diff(Limb.(sess).RWristY(:,2))] ];


end


beh_all.vel_rwrist = [ beh_all.vel_rwristX(:,1), sqrt( beh_all.vel_rwristX(:,2).^2+beh_all.vel_rwristY(:,2).^2 ) ];



%%
interSess = 1000;
maxt=0;
tm = [];

for jj=1:nSess
sess = session_types{use_sessions(jj)};
totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;

tm = [ tm; Vidtime.(sess).bcam(:,2)+maxt-totalt];

end