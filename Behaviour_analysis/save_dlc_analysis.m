%% transfer data

wt = allData.behaviour.whiskerMI(:,1);
wt = [wt(1)-wt(2)+wt(1); wt];
nt = 1:length(wt);

wamp = []; wang = []; waf = []; wph = []; wsp = [];
for sn = use_sessions
sess = session_types{sn};
wamp = [wamp; dlc.(sess).whisk_amp];
wang = [wang; dlc.(sess).whisk_angle];
waf = [waf; dlc.(sess).whisk_angle_filt];
wph = [wph; dlc.(sess).whisk_phase];
wsp = [wsp; dlc.(sess).whisk_set_point];
end

allData.behaviour.whisker_amplitude = [wt, smooth(wamp(nt),30)];
allData.behaviour.whisker_angle = [wt, wang(nt,:)];
allData.behaviour.whisker_angle_filt = [wt, waf(nt)];
allData.behaviour.whisker_phase = [wt, wph(nt)];
allData.behaviour.whisker_setpoint = [wt, wsp(nt)];

%% Linear Regression


allX = { abs( allData.behaviour.speed), allData.behaviour.speed, allData.behaviour.whiskerMI, allData.behaviour.whisker_amplitude, allData.behaviour.whisker_angle_filt, allData.behaviour.whisker_setpoint};
pc1 = allAnalysed.PCA.dff.all.proj(1,:)/allAnalysed.PCA.dff.all.eig_val(1);
pc2 = allAnalysed.PCA.dff.all.proj(2,:)/allAnalysed.PCA.dff.all.eig_val(2);
allY = [allData.neurons.f, pc1', pc2'];
yt = [allData.neurons.time, nanmean(allData.neurons.time,2), nanmean(allData.neurons.time,2)];

nNeu = size(allY,2);
fs = 10; nCV = 10;
lm = [0,0.1, 0.2, 0.5, 1, 5, 10];
evRes = nan( nNeu, length(lm) );
evRes2=evRes;

fld='lincomb'
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld = 'wmi'
allX2 = {allData.behaviour.whiskerMI};
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld='speed'
allX2 = {allData.behaviour.speed};
allX2{1}(:,2)=smooth(allX2{1}(:,2),50);
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld='loco'
allX2 = {abs(allData.behaviour.speed)};
allX2{1}(:,2)=smooth(allX2{1}(:,2),50);
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld='wangle'
allX2 = {allData.behaviour.whisker_angle_filt};
allX2{1}(:,2)=smooth(allX2{1}(:,2),30);
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld='wamp'
allX2 = {allData.behaviour.whisker_amplitude};
allX2{1}(:,2)=smooth(allX2{1}(:,2),30);
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);

fld='wsp'
allX2 = {allData.behaviour.whisker_setpoint};
for jj=1:length(lm)
lambda2 = lm(jj);
[ev, ev2] = get_cv_linear_fit( allY, yt, allX2, fs,  lambda2, nCV );
evRes(:,jj) = nanmean( ev,2);
evRes2(:,jj) = nanmean( ev2,2);
end
allAnalysed.LinModel.neurons.(fld).train = evRes2(1:nNeu-2,:);
allAnalysed.LinModel.neurons.(fld).test  = evRes(1:nNeu-2,:);
allAnalysed.LinModel.pcs.(fld).train = evRes2(nNeu-1:end,:);
allAnalysed.LinModel.pcs.(fld).test  = evRes(nNeu-1:end,:);
