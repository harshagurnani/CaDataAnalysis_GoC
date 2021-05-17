%plot correlations
sess = session_types{use_sessions(1)};
C = Corr.pw.(sess).all.corr(:,:,41);
C(eye(size(C,1))==1)=1;
[B,twom] = modularity(C);
for jj=1:10,[S(1:p.nROIs,jj),Q]= iterated_genlouvain(B); end
S = median(S,2);
[~, cid] = sort(S,'descend');

figure; colormap(redblue)
epoch_type = {'all','loco','whisk'};
for jj=1:nSess
for kk = 1:length(epoch_type)
    ep = epoch_type{kk};
    sess = session_types{use_sessions(jj)};
    subplot(3,nSess,(kk-1)*nSess + jj)
    try    
        nr = floor(size(Corr.pw.(sess).(ep).corr,3)/2);
    imagesc( Corr.pw.(sess).all.corr(cid,cid,nr+1), [-1 1] ); 
    catch
    end
    title([strrep(sess,'_',' '), ' - ', ep])
end
end

nrows = nSess;
figure; 
for jj=1:nSess
    sess = session_types{use_sessions(jj)};
    subplot(2,nrows,jj)
    try    
    nr = floor(size(Corr.pw.(sess).all.corr,3)/2);
    all_pw = Corr.pw.(sess).all.corr(:,:,nr+1);
    all_sig = Corr.pw.(sess).all.corr_sig(:,:,nr+1);
    all_w = Corr.pw.(sess).whisk.corr;
    w_sig = Corr.pw.(sess).whisk.corr_sig;
    id = [abs(all_sig(:))==1 & abs(w_sig(:))==1]; 
    scatter( all_pw(id), all_w(id),'filled' ); 
%     scatter( all_pw(:), all_w(:),'filled' ); 

    catch
    end
    hold on;
    plot([-.3 1],[-.3 1],'k--')
    scatter(nanmedian(all_pw(id)), nanmedian(all_w(id)),'ko','filled')
%     scatter(nanmedian(all_pw(:)), nanmedian(all_w(:)),'ko','filled')
    title([strrep(sess,'_',' '), ' - ', ep])
    subplot(2, nrows, nrows+jj)
    nb=5;bins = [-0.05*nb:0.05:0.05*nb];
    hist( all_w(id) - all_pw(id), bins)
    [c,b]=hist( all_w(id) - all_pw(id), bins);
    hold on;
    plot([0 0 ], [0,1.2*c(nb+1)],'k--')
    scatter(nanmedian(all_w(id) - all_pw(id)), 1.1*c(nb+1), 'o','filled')
    xlim([-.05*nb,0.05*nb])
end

%% all pairwise - combined
sess = 'combined';
C = Corr.pw.(sess).all.corr(:,:,41);
C(eye(size(C,1))==1)=1;
[B,twom] = modularity(C);
for jj=1:10,[S(1:p.nROIs,jj),Q]= iterated_genlouvain(B); end
S = median(S,2);
[~, cid] = sort(S,'descend');

figure;  colormap(redblue);
subplot(2,1,1)
nr = floor(size(Corr.pw.(sess).all.corr,3)/2);
pw=Corr.pw.(sess).all.corr(cid,cid,nr+1);
pw(eye(p.nROIs)==1)=1;
imagesc( pw, [-1 1] ); 
hp4 = get(subplot(2,1,1),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.02  hp4(2)/1.75]);


subplot(2,1,2)
id = (abs(Corr.pw.(sess).all.corr_sig(cid,cid,nr+1))==1);
[c,b]=hist(pw(id),[-1:0.1:1]);
a=area(b,c);
a.FaceAlpha = 0.5;
hold on; plot([0 0],[0 1.2*max(c)],'k--')
plot([median(pw(id)), median(pw(id))], [0 max(c)], 'b','linewidth',1.5)



%% combined corr
figure;
sess = 'combined';
C = Corr.pw.(sess).all.corr(:,:,41);
C(eye(size(C,1))==1)=1;
[B,twom] = modularity(C);
for jj=1:10,[S(1:p.nROIs,jj),Q]= iterated_genlouvain(B); end
S = median(S,2);
[~, cid] = sort(S,'descend');


colormap(redblue)
epoch_type = {'all','whisk'};
ttl = {'r_0','r_w'};
for kk = 1:length(epoch_type)
    ep = epoch_type{kk};
    subplot(2,2,kk)
    try    
    nr = floor(size(Corr.pw.(sess).(ep).corr,3)/2);
    imagesc( Corr.pw.(sess).all.corr(cid,cid,nr+1), [-1 1] ); 
    catch
    end
    title([ep,' (', ttl{kk},')'])
end
hp4 = get(subplot(2,2,2),'Position')
colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)/1.9])


subplot(2,2,3)
try    
nr = floor(size(Corr.pw.(sess).all.corr,3)/2);
all_pw = Corr.pw.(sess).all.corr(:,:,nr+1);
all_sig = Corr.pw.(sess).all.corr_sig(:,:,nr+1);
all_w = Corr.pw.(sess).whisk.corr;
w_sig = Corr.pw.(sess).whisk.corr_sig;
id = [abs(all_sig(:))==1 & abs(w_sig(:))==1]; 
scatter( all_pw(id), all_w(id),'filled' ); 
%     scatter( all_pw(:), all_w(:),'filled' ); 

catch
end
hold on;
plot([-.2 1],[-.2 1],'k--')
xlim([-.2 1])
ylim([-.2 1])
scatter(nanmean(all_pw(id)), nanmean(all_w(id)),'ko','filled')
%     scatter(nanmedian(all_pw(:)), nanmedian(all_w(:)),'ko','filled')
title([strrep(sess,'_',' '), ' - ', ep])
ylabel('Corr during whisking r_w')
xlabel('All corr r_0')

subplot(2,2, 4)
nb=7;bins = [-0.05*nb:0.05:0.05*nb];
[c,b]=hist( all_w(id) - all_pw(id), bins);
a=area( b, c/sum(c));
a.FaceAlpha=0.3;
hold on;
plot([0 0 ], [0,max(c)/sum(c)],'k--')
scatter(nanmedian(all_w(id) - all_pw(id)), 1.1*c(nb+1)/sum(c), 'ko','filled')
xlim([-.05*nb,0.05*nb])
xlabel('r_w - r_0')
ylabel('Frac of cell pairs')

%%
nComp = size(activity_all.pc.proj,1);



figure; nComp=size(activity_all.pc.proj,1);
for jj=1:size(activity_all.pc.proj,1)
   plot(  activity_all.t(:,1)/1000, ((prctile(activity_all.pc.eigvec(:,jj),20)>-0.15)*2-1)*activity_all.pc.proj(jj,:) + 10*(nComp-jj+1)); 
   hold on
    
end

maxt=0;


ispuff=false;
% Plot top Pc with behaviour
nComp = size(activity_all.pc.proj,1);
figure; hold on;cj=colormap(hot(2*nComp));
cj=cj(3:end,:);
nComp_plot=min(nComp,  5);
ispuff=false;
nPuff = size(concat_all_events.AP);
for trial=1:nPuff
    ps = concat_all_events.AP(trial,1)/1000; pe = concat_all_events.AP(trial,2)/1000; %second
    a = area( [ps,pe], [10*nComp+20 10*nComp+20], -10)   ;
    a.FaceColor = 'm'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;a.EdgeColor = 'r';
    beh(1)=a;
    ispuff=true;
end
nPuff = size(concat_all_events.puff_out);
for trial=1:nPuff
    ps = concat_all_events.puff_out(trial,1)/1000; pe = concat_all_events.puff_out(trial,2)/1000; %second
    a = area( [ps,pe], [10*nComp+20 10*nComp+20], -10)   ;
    a.FaceColor = 'b'; a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;a.EdgeColor = 'k';
    beh(1)=a;
    ispuff=true;
end
     
for jj=1:nComp_plot
   plot(  activity_all.t(:,1)/1000, ((prctile(activity_all.pc.eigvec(:,jj),10)>-0.2)*2-1)*activity_all.pc.proj(jj,:) + 10*(nComp_plot-jj),'color', cj(jj,:)); 
   hold on
    
end

beh(2)=plot(beh_all.loco(:,1)/1000, beh_all.loco(:,2)+10*nComp_plot+20, 'b','linewidth', 1)
beh(3)=plot(beh_all.whisk(:,2)/1000, 20*beh_all.whisk(:,1)+25+10*nComp_plot, 'g','linewidth', 1)
if ispuff,legend(beh,'Air Puff','Wheel speed', 'Whisking MI')
else
    legend(beh(2:3),'Wheel speed', 'Whisking MI')
end

% plot soma weights
use_pc = [1 2 3];
soma_wts = activity_all.pc.eigvec(:,use_pc);
figure;
scatter3( soma_wts(:,1), soma_wts(:,2), soma_wts(:,3), 'filled')
xlabel(['PC', num2str(use_pc(1))])
ylabel(['PC',num2str(use_pc(2))])
zlabel(['PC', num2str(use_pc(3))])

%normalised soma wts
use_pc = [2 3 4];
soma_wts = activity_all.pc.eigvec(:,use_pc);
figure;
scatter3( soma_wts(:,1), soma_wts(:,2)./soma_wts(:,1), soma_wts(:,3)./soma_wts(:,1), 'filled')
xlabel(['PC', num2str(use_pc(1))])
ylabel(['PC',num2str(use_pc(2))])
zlabel(['PC', num2str(use_pc(3))])
ylim([-5 5])
zlim([-5 5])