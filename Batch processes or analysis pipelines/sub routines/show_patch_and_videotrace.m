%% Plotting parameters
[~,zid] = sort(p.POI_coordinates.ROIs_Z(:,1)); zid = ([zid(1:2:p.nROIs); zid(2:2:p.nROIs)]);  % sort cells by depth
nr = ceil(p.nROIs/2); nr2 = p.nROIs-nr;


% Axis settings
maxt = 10*ceil(max(Norm.(sess).t(:)/1000)/10)+10;                         % Imaging length (time in ms)
startt = 10;                                                                % number of dff timepoints to not plot (because of causal smooth - omit first startt points)
dispt = 2;                                                                  % time in sec % how much i'm shifting ALL traces by to get whitespace.
t1 = dispt; t2 = maxt+dispt;                                                % starting x (seconds) of two columns of traces

scale_refx = 10;                                                            % Scaled width of patch mean images  (pixels -> seconds on x-axis)
scale_refy = 50;                                                            % Scaled height of patch mean images (pixels -> scaled dff on y-axis)
scale_dff = 30;                                                             % Scale ( for dF/F )
graylim = [150 800];                                                        % contrast scale for mean patch image

trace_y = scale_refy*nr;
speed_y = scale_refy*nr + 10;
whisk_y = scale_refy*(nr+1)+10;
min_y = -10;


% color settings



%for traces
if ~exist( 'cj', 'var')
cjet = colormap( jet(p.nROIs) ); close all;
if nr2>0, cj = [cjet(1:2:p.nROIs,:); cjet(2:2:p.nROIs,:)];                        % Colormap for traces sorted by depth
else,     cj = cjet;
end
clear cjet
end

% behaviour colours
loco_color = 'b';  % locomotion
whisk_color= 'g';  % whisking
clear leg
%% Plotting

figure;
%for video
colours = colormap(jet(20));

%--------       P1 -  Plot mean image patches
colormap gray
hold on; for jj=1:nr , roi=zid(jj); imagesc([-scale_refx,0], [scale_refy*(jj-1),scale_refy*(jj-1)], Seg.mean_image{roi},graylim); end
hold on; for jj=1:nr2, roi=zid(nr+jj); imagesc([maxt-scale_refx,maxt], [scale_refy*(jj-1),scale_refy*(jj-1)], Seg.mean_image{roi},graylim); end

%-------       P2 - Plot behaviour traces
if get_beh_epoch
if isfield( Loco, 'Detected')
    if isfield(Loco.Detected, sess)
    % plot running modules 
    nRuns = size(Loco.Detected.(sess), 1); 
    for pp = 1:nRuns
       a=area(t1+[Loco.Detected.(sess)(pp,1), Loco.Detected.(sess)(pp,2)]/1000, [speed_y, speed_y], min_y);
       a.FaceColor = loco_color; a.FaceAlpha = 0.2; a.EdgeAlpha = 0.2; a.EdgeColor = loco_color;
       a = area( t2+[Loco.Detected.(sess)(pp,1), Loco.Detected.(sess)(pp,2)]/1000, [speed_y,speed_y], min_y)   ;
       a.FaceColor = loco_color; a.FaceAlpha = 0.2; a.EdgeAlpha = 0.2; a.EdgeColor = loco_color;
    end
    end
end

if isfield( Whisk, 'Detected')
    % plot whisking modules 
    nRuns = size(Whisk.Detected.(sess), 1); 
    for pp = 1:nRuns
       a=area(t1+[Whisk.Detected.(sess)(pp,1), Whisk.Detected.(sess)(pp,2)]/1000, [whisk_y, whisk_y], min_y );
       a.FaceColor = whisk_color; a.FaceAlpha = 0.2; a.EdgeAlpha = 0.2;a.EdgeColor = whisk_color;
       a = area( t2+[Whisk.Detected.(sess)(pp,1), Whisk.Detected.(sess)(pp,2)]/1000, [whisk_y, whisk_y], min_y )   ;
       a.FaceColor = whisk_color; a.FaceAlpha = 0.2; a.EdgeAlpha = 0.2;a.EdgeColor = whisk_color;
    end
end
end

% Add airpuff if needed
if add_puff
    for puffnum = 1:length(Puff.(sess))  %trial=1:p.nTrials
        %ps = pufftime.(sess)(trial)/1000; pe = (pufftime.(sess)(trial)+Puff.dur)/1000; %second
        ps = Puff.(sess)(puffnum,1)/1e3; pe = Puff.(sess)(puffnum,2)/1e3;
        
        a = area( [dispt+ps,dispt+pe], [50*nr,50*nr], -10)   ;
        a.FaceColor = 'm'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;a.EdgeColor = 'r';
        a = area( [maxt+dispt+ps,maxt+dispt+pe], [50*nr,50*nr], -10)   ;
        a.FaceColor = 'm'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;a.EdgeColor = 'r';
    end
end

% plot dF/F
for jj=1:nr, roi=zid(jj);plot(Norm.(sess).t(startt:end,roi)/1000+dispt, scale_dff*(Norm.(sess).n(startt:end,roi)-min(Norm.(sess).n(startt:end,roi)))+(jj-1)*scale_refy+10,'color',cj(jj,:),'linewidth',1); end
for jj=1:nr2, roi=zid(nr+jj);plot(Norm.(sess).t(startt:end,roi)/1000+maxt+dispt, scale_dff*(Norm.(sess).n(startt:end,roi)-min(Norm.(sess).n(startt:end,roi)))+(jj-1)*scale_refy+10,'color',cj(jj+nr,:),'linewidth',1); end
xlim([-15 2*maxt+10])
ylim([2*min_y scale_refy*(nr+1)+scale_refy*nVid+scale_refy])
%Plot scale bars for 100%dFF and 10s
plot( [2*maxt,2*maxt], [scale_refy*(nr-1), scale_refy*(nr-1)+scale_dff], 'k', 'LineWidth', 1.5)
plot( [2*maxt-10,2*maxt], [scale_refy*(nr-1), scale_refy*(nr-1)], 'k', 'LineWidth', 1.5)


% Plot speed
if get_speed
    if ~isempty(Speed.(sess))
    leg(1)=plot(Speed.(sess)(:,1)/1000+dispt, 3*smooth(abs(Speed.(sess)(:,2)),101)+50*nr+50,'k', 'Linewidth',1.5);
    plot(Speed.(sess)(:,1)/1000+maxt+dispt, 3*smooth(abs(Speed.(sess)(:,2)),101)+50*nr+50,'k', 'Linewidth',1.5)
    else
    leg(1) = plot(NaN, NaN, 'k');
    end
else
    leg(1) = plot(NaN, NaN, 'k');    
end

%Plot motion index
if get_videoMI
for vn =1:nVid
   video = exp.Video.Segments{vn};
   if isfield(VidMI.(sess), video)
   leg(1+vn) = plot(VidMI.(sess).(video)(:,2)/1000+dispt, 50*VidMI.(sess).(video)(:,1) +50*nr+50*vn+50, 'Color',colours(vn,:),'Linewidth',1);
   plot(VidMI.(sess).(video)(:,2)/1000+maxt+dispt, 50*VidMI.(sess).(video)(:,1) +50*nr+50*vn+50, 'Color',colours(vn,:),'Linewidth',1)
   end
end
newlegs = cellfun(@(s) strrep(s,'_',' '),exp.Video.Segments,'UniformOutput',false);
legend( leg, 'Speed', newlegs{:})
end

dx=10*floor(maxt/100);
set(gca,'xtick', [t1:dx:maxt-2*dx/3, t2:dx:2*maxt+2-2*dx/3], 'xticklabel',num2str([0:dx:maxt-2*dx/3,0:dx:maxt-2*dx/3]'))
set(gca,'ytick',[])

title( [strrep(sess,'_',' ') ' experiment'])
clear nr nr2 roi jj startt graylim scale