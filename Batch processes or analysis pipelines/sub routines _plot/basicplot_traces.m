function [plotfig, PlotParams] = basicplot_traces ( PlotParams, NormF, NormT, p, sessTitle, Speed, PuffTimes, meanImage )
%% Basic Plot - Traces with Mean Patch Image/ROI

    PlotParams = basicplot_params( PlotParams );
    
    % Get essential Params
    nROI = size(NormF, 2);
    maxt = 10*ceil(max(NormT(:)/1000)/10)+2*PlotParams.patchSize(1);           % second
    

    t1 = PlotParams.TmOffset; t2 = maxt+PlotParams.TmOffset;                                            % starting x (seconds) of two columns of traces

    
    if ~isempty(PlotParams.sortDepth)
    % sort cells by depth
        [~,zid] = sort(p.POI_coordinates.ROIs_Z(:,1), PlotParams.sortDepth); 
    else, zid =1:nROI;
    end
    
    
    % Two Column?
    plotfig = figure;
    if isempty(PlotParams.two_col)
        nr = ceil(nROI/2); 
        if nr< 7, nr  = nROI; PlotParams.two_col = false;   %Plot only one column if too few cells
        else,     PlotParams.two_col = true; end
    end
    if PlotParams.two_col
        nr = ceil(nROI/2);  
        zid = ([zid(1:2:nROI); zid(2:2:nROI)]);
    end
    nr2 = nROI-nr;

    % Colourmap
    if ~isempty( PlotParams.ColorMap )
        cjet = colormap( jet(nROI) ); 
    elseif ~exist( 'cj', 'var') 
        cjet = colormap( jet(nROI) );
    end
    
    % swap colours if two column
    if nr2>0, cj = [cjet(1:2:nROI,:); cjet(2:2:nROI,:)];
    else,     cj = cjet;
    end
    clear cjet


    yjump = PlotParams.patchSize(2);
    
    %% PLOT1: plot mean patch images
    if ~isempty(meanImage)
    %Plot Patches
        colormap gray; hold on;
        for jj=1:nr
        roi=zid(jj);     imagesc( [-PlotParams.patchSize(1),        0],  [PlotParams.patchSize(2)*(jj-1),PlotParams.patchSize(2)*(jj-1)],  meanImage{roi},      PlotParams.graylim );    end
        for jj=1:nr2, 
        roi=zid(nr+jj);  imagesc( [maxt-PlotParams.patchSize(1), maxt],  [PlotParams.patchSize(2)*(jj-1),PlotParams.patchSize(2)*(jj-1)],  meanImage{roi},      PlotParams.graylim );    end
    end
    
    
    %% PLOT2: Add airpuff if needed
    if PlotParams.addPuff
       if ~isempty(PuffTimes)
        for trial=1:nTr
            ps = PuffTimes(trial,1)/1000; pe = PuffTimes(trial,2)/1000;     %second
            a = area( [t1+ps,t1+pe], [PlotParams.patchSize(2),PlotParams.patchSize(2)], -10)   ;
            a.FaceColor = PlotParams.EventColor.Puff;    a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;
            if nr2>0
            a = area( [t2+ps,t2+pe], [PlotParams.patchSize(2),PlotParams.patchSize(2)], -10)   ;
            a.FaceColor = PlotParams.EventColor.Puff;    a.FaceAlpha = 0.3; a.EdgeAlpha = 0.3;
            end
        end
       end
    end

    %% PLOT3: Add dff traces
    for jj=1:nr,  roi=zid(jj);   plot(NormT(PlotParams.TmOffset:end,roi)/1000  +t1,    PlotParams.dff_scale*(NormF(PlotParams.TmOffset:end,roi) -min(NormF(PlotParams.TmOffset:end,roi)))+(jj-1)*yjump+10,'color',cj(jj,:),   'linewidth',1); end
    for jj=1:nr2, roi=zid(nr+jj);plot(NormT(PlotParams.TmOffset:end,roi)/1000  +t2,    PlotParams.dff_scale*(NormF(PlotParams.TmOffset:end,roi) -min(NormF(PlotParams.TmOffset:end,roi)))+(jj-1)*yjump+10,'color',cj(jj+nr,:),'linewidth',1); end


    %% PLOT4: Plot running speed 
    if PlotParams.addSpeed && ~isempty(Speed)
        plot(Speed(:,1)/1000+t1, 3*smooth(abs(Speed(:,2)),101) +PlotParams.patchSize(2)*nr+10, 'k', 'Linewidth',1.5)
        plot(Speed(:,1)/1000+t2, 3*smooth(abs(Speed(:,2)),101) +PlotParams.patchSize(2)*nr+10, 'k', 'Linewidth',1.5)
    end

    % axis limits
    xlim([-15 maxt+10 + maxt*(nr2>0)])
    ylim([-30 PlotParams.patchSize(2)*nr+30])

    % ticks
    dx=10*floor(maxt/100);
    if nr2> 0
        set(gca,'xtick', [PlotParams.TmOffset:dx:maxt-2*dx/3, maxt+PlotParams.TmOffset:dx:2*maxt+2-2*dx/3], 'xticklabel',num2str([0:dx:maxt-2*dx/3,0:dx:maxt-2*dx/3]'))
    else
        set(gca,'xtick', [PlotParams.TmOffset:dx:maxt-2*dx/3], 'xticklabel',num2str([0:dx:maxt-2*dx/3]'))
    end
    set(gca,'ytick',[])

    % scalebar for dff
    if nr2>0, xx = 2*maxt; else, xx = maxt+20; end
    plot( [xx,xx],      [yjump*(nr-1), yjump*(nr-1)+PlotParams.dff_scale], 'k', 'LineWidth', 1.5)
    plot( [xx-10,xx],   [yjump*(nr-1), yjump*(nr-1)], 'k', 'LineWidth', 1.5)

    title( [strrep(sessTitle, '_', ' ') ' experiment'])
    clear nr nr2 roi jj 

end