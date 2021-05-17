if isempty(exp.Use.Masks)
    %% Create masks
    
    if strcmpi(Basic.ROI_type, 'patch'),    Seg.ReThresh = 0.5;
    elseif strcmpi(Basic.ROI_type, 'plane'),Seg.ReThresh = 0.7;
    end    
    Seg.method = 'sobel';
    
    % Find edges and fill holes - to get contrasted convex objects
    for roi=1:nROI
        [~,threshold]=edge(im{roi}, Seg.method);    Seg.threshold(roi) = Seg.ReThresh*threshold;
        BW      = edge(im{roi}, Seg.method, Seg.threshold(roi) );           % re-threshold
        BWdil   = imdilate(BW, strel('square',3));                          % blur/expand
        Seg.fullmasks{roi} = imfill(BWdil, 'holes');%.*keeppix{roi}; Do this multiplication afterwards  % fill in holes
    
        % Make background mask    
        Seg.all_background_masks{roi} =  ~(imdilate(Seg.fullmasks{roi}, strel('disk',4,4))); % inverse of detected ROI components
        Seg.all_background_masks{roi} = Seg.all_background_masks{roi} .*keeppix{roi};

        % Dark background
        bg_lim = prctile(Seg.mean_image{roi}(Seg.fullmasks{roi}),10);       % lowest 10th of fullmask
        Seg.background_masks{roi} = Seg.all_background_masks{roi} .*(Seg.mean_image{roi}<bg_lim);

    end
        
        
    if strcmpi(Basic.ROI_type, 'patch')       
        for roi=1:nROI
            % keep maximal connected components for each patch
            Seg.masks{roi} = false(size(Seg.fullmasks{roi}));
            Seg.masks{roi}( get_bigg_conncomp(Seg.fullmasks{roi}) ) = true;
            Seg.masks{roi} = Seg.masks{roi} .*keeppix{roi};
            Seg.mask_roi(roi,1) = roi;
        end   
        clear BW BWdil threshold

    elseif strcmpi(Basic.ROI_type, 'plane')
        % keep all components greater than given size - in a large patch
        % ('plane')
        if isempty(Basic.Mask_MinSize), Basic.Mask_MinSize = 0; end
        nComp = 0; Seg.masks={};Seg.mask_roi = [];
        for roi=1:nROI
            temp = all_comp_masks(Seg.fullmasks{roi}, Basic.Mask_MinSize );
            nComp = length(Seg.masks);
            Seg.masks(nComp+1:nComp+length(temp)) = temp;
            Seg.mask_roi(nComp+1:nComp+length(temp),1) = deal(roi);
        end
    end
    
    

else
    %% Use masks
    load( exp.Folder.UseMasks, 'Seg' );
    fullmasks = Seg.fullmasks; masks = Seg.masks; mask_roi = Seg.mask_roi;
    background = Seg.background_masks;
    clear Seg
    Seg = struct( 'fullmasks', fullmasks, 'masks', masks, 'background_masks', background, 'mask_roi', mask_roi );
    clear masks fullmasks background mask_roi
end


%% HELPER
%------------------------------------------------------------------------%
function [ masks ] = all_comp_masks( bin_img, min_mask_size )
components = bwconncomp(bin_img);
numPixels      = cellfun(@numel,components.PixelIdxList);
keepComp = find(numPixels>min_mask_size);

for jj=1:length(keepComp)
    masks{jj} = false(size(bin_img));
    masks{jj}( components.PixelIdxList{keepComp(jj)} );
end
end

function pix = get_bigg_conncomp( bin_img )
components     = bwconncomp(bin_img);
numPixels      = cellfun(@numel,components.PixelIdxList);
[~,idx]        = max(numPixels);
pix            = components.PixelIdxList{idx(1)};
end
