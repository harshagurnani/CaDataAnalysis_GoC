gap = 1;    %second
c=colormap(jet(p.nROIs));
for sn=1:nSess
    sess = session_types{use_sessions(sn)};
   
    if ~isempty(Puff.(sess))
    
    max_trialT = PuffResponses.(sess).AlignedT(end)/1000; %second
    min_trialT = PuffResponses.(sess).AlignedT(1)/1000; %second

    %% Figure 1 - all trial responses
    
    figure;
    % air puff - bars
    for trial=1:p.nTrials
        hold on; 
        b=bar([0 Puff.dur/1000] +(trial-1)*(max_trialT+gap), [0 p.nROIs+5] ); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2; 
    end
    
    % traces
    for trial=1:p.nTrials    
    for roi=1:p.nROIs
        plot( PuffResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), 1.5*PuffResponses.(sess).AlignedN(:, trial, roi) +roi,    'linewidth',1.5,'color',c(roi,:)); 
    end
    end
    
    % speed data
    for trial=1:p.nTrials
        plot( PuffResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), -0.2*PuffResponses.(sess).AlignedLoco(:,trial) + p.nROIs+2,   'b','linewidth',1.5); 
    end
    
    % whisking data
    for trial=1:p.nTrials
        plot( PuffResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), PuffResponses.(sess).AlignedWhisk(:,trial) + p.nROIs+3,   'g','linewidth',1.5); 
    end
    
    % trial averaged response
    jj=2+p.nTrials;
    for roi=1:p.nROIs
        plot( PuffResponses.(sess).AlignedT/1000 + jj*(max_trialT+gap), 1.5*nanmean(PuffResponses.(sess).AlignedN(:,:,roi),2)+roi,     'linewidth',1.5,'color',c(roi,:)); 
    end
    b=bar([0 Puff.dur/1000] + jj*(max_trialT+gap), [0 p.nROIs+10] ); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2; 
    ylim([0 p.nROIs+5])
    xlim([min_trialT - gap, (jj+1)*(max_trialT+gap) ])
    ylabel('Cell #'); 
    title(['All AP responses during session ', strrep(sess,'_',' ')])
    savefig([exp.analysed,'\',sess,'_all_AP_responses.fig'])
    
    
    %% Trial averaged response
    xt=[];for roi=1:p.nROIs, xt = [xt,nanmean(PuffResponses.(sess).AlignedN(:,:,roi),2)]; end %trial-averaged dff
    figure;
    imagesc(zscore(xt)');             % plot z-scored per roi
    maxt=size(xt,1);hold on
    
    % Axis tick labels
    % Model: y = a*x + b
    % First point: 
    y1=PuffResponses.(sess).AlignedT(1); x1 = 1;
    % Last point:  
    y2 =PuffResponses.(sess).AlignedT(end); x2 = maxt;
    
    a = (y1-y2)/(x1-x2); b = y1-a*x1;
    x0 = -b/a; x100 = (100-b)/a;
    
    plot([x0 x0],[0 p.nROIs+1], 'k')
    plot([x100 x100],[0 p.nROIs+1], 'k')

    xticky = [-400:200:1500];
    xtickx = (xticky-b)/a;
    
    set(gca,'xtick',xtickx,'xticklabels',num2str(xticky')  )
    ylabel('Cell #'); xlabel('Time (ms)')
    title(['Average response to AP during session ', strrep(sess,'_',' ')])
    savefig([exp.analysed,'\',sess,'_avg_AP_response.fig'])

    end 
end