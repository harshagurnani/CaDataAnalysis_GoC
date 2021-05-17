interSess = 1000;
maxt=0;


%%
nVids = numel(exp.Video.Segments);
for jj=1:nVids
   video = exp.Video.Segments{jj};
   beh_all.(video) = [];
end

for jj=1:nSess
sess = session_types{use_sessions(jj)};
totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;

for kk=1:nVids
   video = exp.Video.Segments{kk};
   [~,id]=unique(VidMI.(sess).(video)(:,2));
   id = sort(id,'ascend');
   beh_all.(video) = [beh_all.(video); VidMI.(sess).(video)(id,:)+[0 maxt-totalt] ];
end
end