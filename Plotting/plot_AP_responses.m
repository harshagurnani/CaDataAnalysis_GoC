for sn=1:nSess
    sess = session_types{use_sessions(sn)};
   
    if ~isempty(Puff.(sess))
    
    figure;
    %% all responses
    gap = 1;c=colormap(jet(p.nROIs));
    for jj=1:p.nTrials,hold on; b=bar([0.5 0.6]+PuffResponses.(sess).t(1,1,1)/1000+(jj-1)*((PuffResponses.(sess).t(end,1,end)-PuffResponses.(sess).t(1,1,end))/1000+gap), [0 40]); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2; end
    for jj=1:p.nTrials
    hold on;
    for roi=1:p.nROIs
    plot(squeeze(PuffResponses.(sess).t(:,1,roi))/1000+(jj-1)*((PuffResponses.(sess).t(end,1,roi)-PuffResponses.(sess).t(1,1,roi))/1000+gap), 1.5*(squeeze(PuffResponses.(sess).n(:,jj,roi))-nanmean(PuffResponses.(sess).n(1:10,jj,roi)))+roi,'linewidth',1.5,'color',c(roi,:)); end
    end
    % speed data
    for jj=1:p.nTrials
    plot(PuffResponses.(sess).locot{1}/1000+(jj-1)*((PuffResponses.(sess).locot{1}(end)-PuffResponses.(sess).locot{1}(1))/1000+gap), p.nROIs+2-0.2*smooth(PuffResponses.(sess).loco{jj},1),'k','linewidth',1.5); end
    for jj=1:p.nTrials
    plot(PuffResponses.(sess).whiskt{1}/1000+(jj-1)*((PuffResponses.(sess).whiskt{1}(end)-PuffResponses.(sess).whiskt{1}(1))/1000+gap), p.nROIs+3+2*smooth(PuffResponses.(sess).whisk{jj},1),'b','linewidth',1.5); end
    for roi=1:p.nROIs
    plot(squeeze(PuffResponses.(sess).t(:,1,roi))/1000+(p.nTrials+1)*((PuffResponses.(sess).t(end,1,roi)-PuffResponses.(sess).t(1,1,roi))/1000+gap), 1.5*nanmean(PuffResponses.(sess).n(:,:,roi)-nanmean(PuffResponses.(sess).n(1:5,:,roi),1),2)+roi,'linewidth',1.5,'color',c(roi,:)); end
    jj=2+p.nTrials;hold on; b=bar([0.5 0.6]+PuffResponses.(sess).t(1,1,1)/1000+(jj-1)*((PuffResponses.(sess).t(end,1,end)-PuffResponses.(sess).t(1,1,end))/1000+gap), [0 40]); b.FaceColor='m';b.FaceAlpha=0.3;b.EdgeAlpha=0.2;
    ylim([0 p.nROIs+6])
    ylabel('Cell #'); 
    title(['All AP responses during session ', strrep(sess,'_',' ')])
    savefig([exp.analysed,'\',sess,'_all_AP_responses.fig'])
    
    
    %% averaged response
    xt=[];for roi=1:p.nROIs, xt = [xt,nanmean(PuffResponses.(sess).n(:,:,roi)-nanmean(PuffResponses.(sess).n(1:5,:,roi),1),2)]; end
    figure;imagesc(zscore(xt)');maxt=size(xt,1);hold on

    % Model: y = a*x + b
    % First point: 
    y1=PuffResponses.(sess).t(1)-Puff.(sess)(1); x1 = 1;
    % Last point:  
    y2 =PuffResponses.(sess).t(maxt,1,1)-Puff.(sess)(1); x2 = maxt;
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