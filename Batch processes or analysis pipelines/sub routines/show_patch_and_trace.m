% sort cells by depth
[~,zid] = sort(p.POI_coordinates.ROIs_Z(:,1)); zid = ([zid(1:2:nROI); zid(2:2:nROI)]);

% plot settings
nr = ceil(nROI/2); 
if nr< 7, nr  = nROI; end   %Plot only one column if too few cells

% nr=nROI;

nr2 = nROI-nr;

if ~exist( 'cj', 'var')
cjet = colormap( hot(floor(2*nROI)) ); close all;
if nr2>0, cj = [cjet(1:2:nROI,:); cjet(2:2:nROI,:)];
else,     cj = cjet;
end
clear cjet
end

maxt = 10*ceil(max(Seg.(sess).t(:)/1000)/10)+10;  % second
scale = 80;                                         % factor to multiply dff 
graylim = [50 800];                                 % mean patch image - range
startt = 10;                                        % how many timepoints to drop from the beginning?
dispt = 2;                                          % second% shift traces along x-axis (space between patch and dff trace)

t1 = dispt; t2 = maxt+dispt;                        % starting x (seconds) of two columns of traces

%% plot mean patch images
figure;
colormap gray
hold on; for jj=1:nr ,
roi=zid(jj);     imagesc([-10,0], [50*(jj-1),50*(jj-1)], Seg.mean_image{roi},        graylim);    end
hold on; for jj=1:nr2, roi=zid(nr+jj);  imagesc([maxt-10,maxt], [50*(jj-1),50*(jj-1)], Seg.mean_image{roi}, graylim);    end

%% Add airpuff if needed
if add_puff
   if ~isempty(Puff.(sess))
    for trial=1:nTr
        ps = pufftime.(sess)(trial)/1000; pe = (pufftime.(sess)(trial)+Puff.dur)/1000; %second
        a = area( [t1+ps,t1+pe], [50*nr,50*nr], -10)   ;
        a.FaceColor = 'm'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;
        if nr2>0
        a = area( [t2+ps,t2+pe], [50*nr,50*nr], -10)   ;
        a.FaceColor = 'm'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;
        end
    end
   end
end

%% Add dff traces
for jj=1:nr, roi=zid(jj);    plot(Norm.(sess).t(startt:end,roi)/1000  +t1,    scale*(Norm.(sess).n(startt:end,roi) -min(Norm.(sess).n(startt:end,roi)))+(jj-1)*50+10,'color',cj(jj,:),'linewidth',1); end
for jj=1:nr2, roi=zid(nr+jj);plot(Norm.(sess).t(startt:end,roi)/1000  +t2,    scale*(Norm.(sess).n(startt:end,roi) -min(Norm.(sess).n(startt:end,roi)))+(jj-1)*50+10,'color',cj(jj+nr,:),'linewidth',1); end


%% Plot behaviour 
if get_speed
    if ~isempty(Speed.(sess))
    plot(Speed.(sess)(:,1)/1000+t1, 3*smooth(abs(Speed.(sess)(:,2)),101)+50*nr+10,'k', 'Linewidth',1.5)
    plot(Speed.(sess)(:,1)/1000+t2, 3*smooth(abs(Speed.(sess)(:,2)),101)+50*nr+10,'k', 'Linewidth',1.5)
    end
end

% axis limits
xlim([-15 maxt+10 + maxt*(nr2>0)])
ylim([-30 50*(nr+1)])

% ticks
dx=10*floor(maxt/100);
if nr2> 0
    set(gca,'xtick', [dispt:dx:maxt-2*dx/3, maxt+dispt:dx:2*maxt+2-2*dx/3], 'xticklabel',num2str([0:dx:maxt-2*dx/3,0:dx:maxt-2*dx/3]'))
else
    set(gca,'xtick', [dispt:dx:maxt-2*dx/3], 'xticklabel',num2str([0:dx:maxt-2*dx/3]'))
end
set(gca,'ytick',[])

% scalebar for dff
if nr2>0, xx = 2*maxt; else, xx = maxt+20; end
plot( [xx,xx],      [50*(nr-1), 50*(nr-1)+scale], 'k', 'LineWidth', 1.5)
plot( [xx-10,xx],   [50*(nr-1), 50*(nr-1)], 'k', 'LineWidth', 1.5)

title( [strrep(sess, '_', ' ') ' experiment'])
clear nr nr2 roi jj startt graylim scale