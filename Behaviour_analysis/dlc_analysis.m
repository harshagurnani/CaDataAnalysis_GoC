%% whisker pad
wpad = wm.Variables
ws = wpad(:,[1,2]);
wl = wpad(:,[3,4]);
vec = wl-ws;
vec = vec./sqrt( sum(vec.*vec,2));
norm(vec(1,:));
vec = mean(vec,1);
vec = vec./sqrt( sum(vec.*vec,2));
wpad_vm = vec;
plot([0,wpad_vm(1)]*500, [0, wpad_vm(2)]*500,'linewidth', 2); hold on;


%% whiskers
nw = 4;
allsess = fieldnames(w);
for sn = 1:length(allsess)
    sess = allsess{sn};
    figure; suptitle(sess)
    w1 = w.(sess).Variables;
    w1(:,1)=[];  w1( :, 3*(1:3*nw) ) = [];
    whiskers = cell(nw,1);
    for ww=1:nw, whiskers{ww} =  w1(:, (ww-1)*6 + (1:6) ); end
    for ww=1:nw,  for jj=1:3, scatter( whiskers{ww}( :, 1 + (jj-1)*2 ), whiskers{ww}( :, 2 + (jj-1)*2 ) ); hold on; end
    end
    plot([0,wpad_vm(1)]*500, [0, wpad_vm(2)]*500,'linewidth', 2); hold on;
    
    whisker_angle
    dlc.(sess).wpad = wpad_vm;
    dlc.(sess).whiskers = whiskers;
    dlc.(sess).whisk_angle = angle;
    dlc.(sess).whisk_skeleton = vm;
    
    
end
plot(zscore(angle)/8+[1:4])

%% variables for 1 whisker
for sn = 1:length(allsess)
    sess = allsess{sn};
    dlc.(sess).whisk_used = 1;
    [dlc.(sess).whisk_angle_filt,dlc.(sess).whisk_set_point,dlc.(sess).whisk_amp,dlc.(sess).whisk_phase] = get_whisking_variables( dlc.(sess).whisk_angle(:,dlc.(sess).whisk_used), 300 );
end


%% time files
path = 'D:\Work\OneDrive - University College London\backupdata\tracking\hg03 exp03a\';
for sn = 1:length(allsess)
    sess = allsess{sn};
    fid = fopen( [path,tf.(sess),'\WhiskersCam-relative times.txt']);
    t = textscan( fid, '%d %d');
    dlc.(sess).time = [double(t{1}), double(t{2})];
end