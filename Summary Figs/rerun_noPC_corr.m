maxt=0;
interSess = 1000;
for jj=1:nSess
sess = session_types{use_sessions(jj)};
totalt = Norm.(sess).t(end-1); maxt = maxt+totalt+interSess;
end
nROI= p.nROIs;
x = activity_all.soma;
xnew = subspace_svd( x, 2:nROI);
xnew2 = subspace_svd( x, [1,3:nROI]);
xnew3 = subspace_svd( x, [1,2,4:nROI]);
[Corr.pw.noPC3.all.lags, Corr.pw.noPC3.all.corr, ~, Corr.pw.noPC3.all.conf_interval, Corr.pw.noPC3.all.corr_sig ] = ...
all_period_lagcorr_pairwise_dff_new( [0 maxt], xnew3, activity_all.t, ...
Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
[Corr.pw.noPC2.all.lags, Corr.pw.noPC2.all.corr, ~, Corr.pw.noPC2.all.conf_interval, Corr.pw.noPC2.all.corr_sig ] = ...
all_period_lagcorr_pairwise_dff_new( [0 maxt], xnew2, activity_all.t, ...
Corr.pw.params.ds_rate, Corr.pw.params.maxlag, [], [], Corr.pw.params.nshuffle, false );
allc=[];
for jj=1:nROI, for kk=jj+1:nROI,
allc=[allc; Corr.pw.combined.all.corr(jj,kk,41)];
end
end
figure;[n0,b]=hist(allc,[-1:0.05:1]);
plot(b,n0)
hold on;
allc=[];
for jj=1:nROI, for kk=jj+1:nROI,
allc=[allc; Corr.pw.noPC1.all.corr(jj,kk,41)];
end
end
[n1,b]=hist(allc,[-1:0.05:1]);
plot(b,n1)
allc=[];
for jj=1:nROI, for kk=jj+1:nROI,
allc=[allc; Corr.pw.noPC2.all.corr(jj,kk,41)];
end
end
[n2,b]=hist(allc,[-1:0.05:1]);
plot(b,n2)
allc=[];
for jj=1:nROI, for kk=jj+1:nROI,
allc=[allc; Corr.pw.noPC3.all.corr(jj,kk,41)];
end
end
[n3,b]=hist(allc,[-1:0.05:1]);
plot(b,n3)