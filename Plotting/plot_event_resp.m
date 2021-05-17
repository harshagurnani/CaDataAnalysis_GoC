function fig = plot_event_resp(EventTimes ,EventResponses, Puff, PCplus)
gap = 1;    %second

sessions = fieldnames(EventResponses);
for sn=1:length(sessions)
    sess = sessions{sn};
    if ~isempty(EventResponses.(sess))
    nTrials = size(EventResponses.(sess).AlignedN,2);
    max_trialT = EventResponses.(sess).AlignedT(end)/1000; %second
    min_trialT = EventResponses.(sess).AlignedT(1)/1000; %second
    nROIs = size(EventResponses.(sess).AlignedN,3);
    if isempty (PCplus), Pcplus = ones(nROIs,1); end
    %% Figure 1 - all trial responses

    fig(sn) = figure;c=colormap(hot(2*nROIs));hold on;
    % air puff - bars
%     if ~isempty(EventTimes.(sess))
%     for ev=1:size(EventTimes.(sess))
%         hold on; 
%         b=bar([EventTimes.(sess)/1000] +(trial-1)*(max_trialT+gap), [0 nROIs+5] ); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2; 
%     end
%     end
    
    % traces
    for trial=1:nTrials    
    for roi=1:nROIs
        plot( EventResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), PCplus(roi)*1.5*zscore(EventResponses.(sess).AlignedN(:, trial, roi)-EventResponses.(sess).AlignedN(1, trial, roi))/10 +roi,    'linewidth',1.5,'color',c(roi,:)); 
    end
    end

    % speed data
    for trial=1:nTrials
        plot( EventResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), -0.2*EventResponses.(sess).loco.AlignedX(:,trial) + nROIs+2,   'b','linewidth',1.5); 
    end

    % whisking data
    for trial=1:nTrials
        plot( EventResponses.(sess).AlignedT/1000 +(trial-1)*(max_trialT+gap), EventResponses.(sess).whisk.AlignedX(:,trial) + nROIs+3,   'g','linewidth',1.5); 
    end

    % trial averaged response
    jj=2+nTrials;
    for roi=1:nROIs
        plot( EventResponses.(sess).AlignedT/1000 + jj*(max_trialT+gap), PCplus(roi)*1.5*zscore(nanmean(EventResponses.(sess).AlignedN(:,:,roi)-EventResponses.(sess).AlignedN(1,:,roi),2))/10+roi,     'linewidth',1.5,'color',c(roi,:)); 
    end
    b=bar([0 Puff.dur/1000] + jj*(max_trialT+gap), [0 nROIs+10] ); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2; 
    ylim([0 nROIs+5])
    xlim([min_trialT - gap, (jj+1)*(max_trialT+gap) ])
    ylabel('Cell #'); 
    title(['All AP responses during session ', strrep(sess,'_',' ')])
%     savefig([exp.analysed,'\',sess,'_all_AP_responses.fig'])


%     %% Trial averaged response
%     xt=[];for roi=1:nROIs, xt = [xt,nanmean(EventResponses.(sess).AlignedN(:,:,roi),2)]; end %trial-averaged dff
%     figure;
%     imagesc(zscore(xt)');             % plot z-scored per roi
%     maxt=size(xt,1);hold on
% 
%     % Axis tick labels
%     % Model: y = a*x + b
%     % First point: 
%     y1=EventResponses.(sess).AlignedT(1); x1 = 1;
%     % Last point:  
%     y2 =EventResponses.(sess).AlignedT(end); x2 = maxt;
% 
%     a = (y1-y2)/(x1-x2); b = y1-a*x1;
%     x0 = -b/a; x100 = (100-b)/a;
% 
%     plot([x0 x0],[0 nROIs+1], 'k')
%     plot([x100 x100],[0 nROIs+1], 'k')
% 
%     xticky = [-400:200:1500];
%     xtickx = (xticky-b)/a;
% 
%     set(gca,'xtick',xtickx,'xticklabels',num2str(xticky')  )
%     ylabel('Cell #'); xlabel('Time (ms)')
%     title(['Average response to AP during session ', strrep(sess,'_',' ')])
%     savefig([exp.analysed,'\',sess,'_avg_AP_response.fig'])


end
end