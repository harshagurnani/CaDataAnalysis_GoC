interSess = 1000;
maxt=0;
beh_all.pupil = [];
beh_all.slow_whiskerMI = [];
%%

for jj=1:nSess
sess = session_types{use_sessions(jj)}
totalt = Norm.(sess).t(end-1); 
maxt = maxt+totalt+interSess;


[~,id]=unique( Pupil.(sess)(:,1) );
id = sort(id,'ascend');

beh_all.pupil = [beh_all.pupil; Pupil.(sess)(id,:) + [maxt-totalt, 0] ];

[~,id]=unique(VidMI.(sess).(Whisk.vidseg)(:,2));
id = sort(id,'ascend');
uwhisk=smooth(VidMI.(sess).(Whisk.vidseg)(:,1),3);
beh_all.slow_whiskerMI = [beh_all.slow_whiskerMI; [uwhisk(id),  VidMI.(sess).(Whisk.vidseg)(id,2)]+ [0 maxt-totalt] ];
end
