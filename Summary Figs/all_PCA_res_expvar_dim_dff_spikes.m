% explained variance by first component

PC.dff.nDim.Spectral = [];
PC.dff.nDim.Peak = [];
PC.dff.nDim.EV80 = [];

PC.dff.ExpVar.Same = [];
PC.dff.ExpVar.CV = [];
PC.dff.ExpVar.Peak = [];

PC.spk.nDim.Spectral = [];
PC.spk.nDim.Peak = [];
PC.spk.nDim.EV80 = [];

PC.spk.ExpVar.Same = [];
PC.spk.ExpVar.CV = [];
PC.spk.ExpVar.Peak = [];

PC.spk.nROI = [];
PC.dff.nROI = [];


for animal=1:nAnimals
for session = 1:nSess
if ~isempty(SummaryData{animal,session})
    if SummaryData{animal,session}.p.nROIs>5%15

    % Explained variance by top two PCs (dFF) - same and cross-validated
    pc1_2=[SummaryData{animal,session}.PCA.all.soma.explained_var(1),...
        SummaryData{animal,session}.PCA.all.soma.explained_var(2)-SummaryData{animal,session}.PCA.all.soma.explained_var(1)];
    PC.dff.ExpVar.Same = [PC.dff.ExpVar.Same;  pc1_2];
    
    res =nanmean(nanmean(SummaryData{animal,session}.PCA.combined.cval.res,1),3);
    pc1_2cv = [res(1), res(2)-res(1)];
    if ~(any(pc1_2cv)<0),     PC.dff.ExpVar.CV = [PC.dff.ExpVar.CV; pc1_2cv];
    else,                     PC.dff.ExpVar.CV = [PC.dff.ExpVar.CV; [NaN NaN]];    
    end
    
    % spectral dimensionality (dff)
    cc = SummaryData{animal,session}.PCA.all.soma.eig_val; %all eigenvalues
    PC.dff.nDim.Spectral = [PC.dff.nDim.Spectral;sum(cc)^2/ sum(cc.^2)];
    % nPC for exp var = 80%
    PC.dff.nDim.EV80 = [PC.dff.nDim.EV80; find(SummaryData{animal,session}.PCA.all.soma.explained_var>0.8,1)];

    % dimensionality as peak of cross validated EV
    aa = nanmean( nanmean( SummaryData{animal,session}.PCA.combined.cval.res, 3),1);
    x = find( aa ==max(aa)); x = x(end);
    PC.dff.nDim.Peak = [PC.dff.nDim.Peak; x ];
    PC.dff.ExpVar.Peak = [PC.dff.ExpVar.Peak; aa(x) ];
    
    PC.dff.nROI = [PC.dff.nROI; length( SummaryData{animal,session}.PCA.all.soma.x_mean) ];
    %-------------------
    
    % Explained variance by top two PCs (spikes) - same and cross-validated
    pc1_2=[SummaryData{animal,session}.PCA.combined.spikes.explained_var(1),...
        SummaryData{animal,session}.PCA.combined.spikes.explained_var(2)-SummaryData{animal,session}.PCA.combined.spikes.explained_var(1)];
    PC.spk.ExpVar.Same = [PC.spk.ExpVar.Same;  pc1_2];
    
    res =nanmean(nanmean(SummaryData{animal,session}.PCA.FR.cval.res,1),3);
    pc1_2cv = [res(1), res(2)-res(1)];
    if ~(any(pc1_2cv)<0),     PC.spk.ExpVar.CV = [PC.spk.ExpVar.CV; pc1_2cv];
    else,                     PC.spk.ExpVar.CV = [PC.spk.ExpVar.CV; [NaN NaN]];    
    end
    
    % spectral dimensionality (spikes)
    cc = SummaryData{animal,session}.PCA.combined.spikes.eig_val; %all eigenvalues
    PC.spk.nDim.Spectral = [PC.spk.nDim.Spectral;sum(cc)^2/ sum(cc.^2)];
    % nPC for exp var = 80%
    PC.spk.nDim.EV80 = [PC.spk.nDim.EV80; find(SummaryData{animal,session}.PCA.combined.spikes.explained_var>0.8,1)];

    % dimensionality as peak of cross validated EV (spikes)
    aa = nanmean( nanmean( SummaryData{animal,session}.PCA.FR.cval.res, 3),1);
    x = find( aa ==max(aa)); x = x(end);
    PC.spk.nDim.Peak = [PC.spk.nDim.Peak; x ];
    PC.spk.ExpVar.Peak = [PC.spk.ExpVar.Peak; aa(x) ];
   
    PC.spk.nROI = [PC.dff.nROI; length( SummaryData{animal,session}.PCA.combined.spikes.x_mean) ];
    
    end
end    
end
end
PC.dff.ExpVar.CV(PC.dff.ExpVar.CV<0)=NaN;
PC.spk.ExpVar.CV(PC.spk.ExpVar.CV<0)=NaN;

%% Plotting

figure; hold on;
for jj=1:length(PC.dff.ExpVar.Same)
   plot([1 2], PC.dff.ExpVar.Same(jj,:) ,'o','markerfacecolor','b','color','b','linewidth',1,'markersize',4) 
   plot([1.3 2.3], PC.dff.ExpVar.CV(jj,:) ,'o','markerfacecolor','k','color','k','linewidth',1,'markersize',4) 
  
   plot([1 1.3], [PC.dff.ExpVar.Same(jj,1) PC.dff.ExpVar.CV(jj,1)], 'k','linewidth',0.3)   
   plot([2 2.3], [PC.dff.ExpVar.Same(jj,2) PC.dff.ExpVar.CV(jj,2)], 'k','linewidth',0.3)
end


bb=bar(1,nanmean(PC.dff.ExpVar.Same(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(2,nanmean(PC.dff.ExpVar.Same(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(1.3,nanmean(PC.dff.ExpVar.CV(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';
bb=bar(2.3,nanmean(PC.dff.ExpVar.CV(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';

clear b
b(1) = bar(NaN, NaN);b(1).FaceAlpha=0.5;b(1).FaceColor=[0 0 1];
b(2) = bar(NaN, NaN);b(2).FaceAlpha=0.5;b(2).FaceColor='k';


errorbar(2,nanmean(PC.dff.ExpVar.Same(:,2)),nanstd(PC.dff.ExpVar.Same(:,2)),'b','marker','.','linewidth',3)
errorbar(1,nanmean(PC.dff.ExpVar.Same(:,1)),nanstd(PC.dff.ExpVar.Same(:,1)),'b','marker','.','linewidth',3)
errorbar(2.3,nanmean(PC.dff.ExpVar.CV(:,2)),nanstd(PC.dff.ExpVar.CV(:,2)),'k','marker','.','linewidth',3)
errorbar(1.3,nanmean(PC.dff.ExpVar.CV(:,1)),nanstd(PC.dff.ExpVar.CV(:,1)),'k','marker','.','linewidth',3)

legend( b, 'Same session', 'Cross-validated')
ylabel('Explained variance')
xlabel('PC')
xticks([1 2])

%% pc1 2 - spikes


figure; hold on;
for jj=1:length(PC.spk.ExpVar.Same)
   plot([1 2], PC.spk.ExpVar.Same(jj,:) ,'o','markerfacecolor','b','color','b','linewidth',1,'markersize',4) 
   plot([1.3 2.3], PC.spk.ExpVar.CV(jj,:) ,'o','markerfacecolor','k','color','k','linewidth',1,'markersize',4) 
  
   plot([1 1.3], [PC.spk.ExpVar.Same(jj,1) PC.spk.ExpVar.CV(jj,1)], 'k','linewidth',0.3)   
   plot([2 2.3], [PC.spk.ExpVar.Same(jj,2) PC.spk.ExpVar.CV(jj,2)], 'k','linewidth',0.3)
end


bb=bar(1,nanmean(PC.spk.ExpVar.Same(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(2,nanmean(PC.spk.ExpVar.Same(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(1.3,nanmean(PC.spk.ExpVar.CV(:,1)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';
bb=bar(2.3,nanmean(PC.spk.ExpVar.CV(:,2)),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';

clear b
b(1) = bar(NaN, NaN);b(1).FaceAlpha=0.5;b(1).FaceColor=[0 0 1];
b(2) = bar(NaN, NaN);b(2).FaceAlpha=0.5;b(2).FaceColor='k';


errorbar(2,nanmean(PC.spk.ExpVar.Same(:,2)),nanstd(PC.spk.ExpVar.Same(:,2)),'b','marker','.','linewidth',3)
errorbar(1,nanmean(PC.spk.ExpVar.Same(:,1)),nanstd(PC.spk.ExpVar.Same(:,1)),'b','marker','.','linewidth',3)
errorbar(2.3,nanmean(PC.spk.ExpVar.CV(:,2)),nanstd(PC.spk.ExpVar.CV(:,2)),'k','marker','.','linewidth',3)
errorbar(1.3,nanmean(PC.spk.ExpVar.CV(:,1)),nanstd(PC.spk.ExpVar.CV(:,1)),'k','marker','.','linewidth',3)

legend( b, 'Same session', 'Cross-validated')
ylabel('Explained variance')
xlabel('PC')
xticks([1 2])
%%
figure;hold on;
for jj=1:length(PC.dff.nDim.Peak),
     plot([1 2 3], [PC.dff.nDim.Peak(jj), PC.dff.nDim.Spectral(jj) PC.dff.nDim.EV80(jj)] ,'o','markerfacecolor','b','color','b','linewidth',1,'markersize',4) 
     plot([1 2 3], [PC.dff.nDim.Peak(jj), PC.dff.nDim.Spectral(jj) PC.dff.nDim.EV80(jj)] , 'k','linewidth',0.3) 
end
bb=bar(1,nanmean(PC.dff.nDim.Peak),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(2,nanmean(PC.dff.nDim.Spectral),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';
bb=bar(3,nanmean(PC.dff.nDim.EV80),0.3);bb.FaceAlpha=0.5;bb.FaceColor='r';
xticks([1 2 3])
xticklabels({'Peak', 'Spectral', 'EV80'})
ylabel('Dimensionality')


%% spikes dimensionality
figure;hold on;
for jj=1:length(PC.spk.nDim.Peak),
     plot([1 2 3], [PC.spk.nDim.Peak(jj), PC.spk.nDim.Spectral(jj) PC.spk.nDim.EV80(jj)] ,'o','markerfacecolor','b','color','b','linewidth',1,'markersize',4) 
     plot([1 2 3], [PC.spk.nDim.Peak(jj), PC.spk.nDim.Spectral(jj) PC.spk.nDim.EV80(jj)] , 'k','linewidth',0.3) 
end
bb=bar(1,nanmean(PC.spk.nDim.Peak),0.3);bb.FaceAlpha=0.5;bb.FaceColor=[0 0 1];
bb=bar(2,nanmean(PC.spk.nDim.Spectral),0.3);bb.FaceAlpha=0.5;bb.FaceColor='k';
bb=bar(3,nanmean(PC.spk.nDim.EV80),0.3);bb.FaceAlpha=0.5;bb.FaceColor='r';
xticks([1 2 3])
xticklabels({'Peak', 'Spectral', 'EV80'})
ylabel('Dimensionality')
title('Spikes')

%% compare dimensionality
figure; hold on;
plot([0,20],[0,20],'k--');
s(1)=scatter( PC.dff.nDim.Peak, PC.spk.nDim.Peak , 30, 'b', 'filled');
s(2)=scatter( PC.dff.nDim.EV80, PC.spk.nDim.EV80 , 30, 'r', 'filled');
s(3)=scatter( PC.dff.nDim.Spectral, PC.spk.nDim.Spectral , 30, 'k', 'filled');
legend(s, 'CV Peak', 'EV 80', 'Spectral')
xlabel('Dimensionality (dff)')
ylabel('Dimensionality (spikes)');


%% normalise dimensionality
figure; hold on;
plot([0,1],[0,1],'k--');
s(1)=scatter( nanmean(PC.dff.nDim.Peak./PC.dff.nROI), nanmean(PC.spk.nDim.Peak./PC.spk.nROI) , 60, 'b','filled');
s(2)=scatter( nanmean(PC.dff.nDim.Spectral./PC.dff.nROI), nanmean(PC.spk.nDim.Spectral./PC.spk.nROI) , 60, 'k','filled');
s(3)=scatter( nanmean(PC.dff.nDim.EV80./PC.dff.nROI), nanmean(PC.spk.nDim.EV80./PC.spk.nROI) , 60, 'r','filled');
scatter( PC.dff.nDim.Peak./PC.dff.nROI, PC.spk.nDim.Peak./PC.spk.nROI , 30, 'b', 'filled', 'MarkerFaceAlpha',0.4);
scatter( PC.dff.nDim.EV80./PC.dff.nROI, PC.spk.nDim.EV80./PC.spk.nROI , 30, 'r', 'filled', 'MarkerFaceAlpha',0.4);
scatter( PC.dff.nDim.Spectral./PC.dff.nROI, PC.spk.nDim.Spectral./PC.spk.nROI , 30, 'k', 'filled', 'MarkerFaceAlpha',0.4);
legend(s, 'CV Peak', 'Spectral', 'EV 80')
xlabel('Normalised Dimensionality (dff)')
ylabel('Normalised Dimensionality (spikes)');


%% normalise dimensionality
figure; hold on;
s2(1)=scatter( PC.dff.nDim.Peak./PC.dff.nROI, PC.dff.ExpVar.Peak , 30, 'b', 'filled');
s2(2)=scatter( PC.spk.nDim.Peak./PC.spk.nROI, PC.spk.ExpVar.Peak , 30, 'm^', 'filled');
legend(s2, 'dff', 'spikes')
xlim([0 0.5])
ylim([0 1])
xlabel('Normalised Dimensionality (CV)')
ylabel('Explained variance (CV)');
