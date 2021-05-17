function plot_patch_in_volume( POI_coordinates, varargin)
    nROI = length(POI_coordinates.ROIs_X);
    zid = 1:nROI;

    % sort by depth?
    if ~isempty(varargin)
        if varargin{1} == 'z'
            [~,zid] = sort(POI_coordinates.ROIs_Z(:,1)); 
        end
       
    end


    figure;
    cjet = colormap( hot(2*nROI) ); 
    hold on;
    
    [ux,uy,uz]=sphere; rx = 5*512/400; rz = 5;ux=ux*rx; uy = uy*rx; uz=uz*rz;
    % spheres
    for jj=1:nROI        
    roi = zid(jj); 
    if roi==26 || roi==36
    s=surf(ux+mean(POI_coordinates.ROIs_X(roi,:)),uy+mean(POI_coordinates.ROIs_Y(roi,:)),uz+mean(POI_coordinates.ROIs_Z(roi,:)));
    
    s.FaceAlpha = 1;s.FaceColor  = cjet(jj,:);  s.EdgeAlpha = 0;
        end
    end
    
    
    % patches
    for jj=1:nROI
    roi = zid(jj);  
    p=fill3([POI_coordinates.ROIs_X(roi,:), POI_coordinates.ROIs_X(roi,[2 1])], ...
           [POI_coordinates.ROIs_Y(roi,1) POI_coordinates.ROIs_Y(roi,1) POI_coordinates.ROIs_Y(roi,2) POI_coordinates.ROIs_Y(roi,2)],...
           repmat(POI_coordinates.ROIs_Z(roi,1),1,4), ones(1,4));
    p.FaceAlpha = 0.2;p.FaceColor  = cjet(jj,:); 
    p.EdgeColor =cjet(jj,:);%'k';%; 
    p.EdgeAlpha =0.2;
    end
    
%     roi numbers
    for jj=1:nROI
       roi = zid(jj);
       if roi==26 || roi==36
       text( mean(POI_coordinates.ROIs_X(roi,:)), mean(POI_coordinates.ROIs_Y(roi,:)), POI_coordinates.ROIs_Z(roi,1), num2str(jj));
        end
    end
end
   