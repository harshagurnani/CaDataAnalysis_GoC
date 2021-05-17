  
%%
% % ALL TP
% a = [];
% for sess=1:exp.nSess
%     sn = session_types{use_sessions(sess)};
%     a = cat(3, a, reshape( Basic.(sn).r, p.nLinePixels,  p.nPatchLines,p.nTimepoints*p.nTrials, p.nROIs));
% end
% % a=reshape( Basic.baseline.r, p.nLinePixels,  p.nPatchLines,p.nTimepoints*p.nTrials, p.nROIs);
% 
% % temporal smoothing
% b = (a(:,:,1:end-4,:)+a(:,:,2:end-3,:)+a(:,:,3:end-2,:)+a(:,:,4:end-1,:)+a(:,:,5:end,:))/5;
% b = (a(:,:,1:end-6,:)+a(:,:,2:end-5,:)+a(:,:,3:end-4,:)+a(:,:,4:end-3,:)+a(:,:,5:end-2,:)+a(:,:,6:end-1,:)+a(:,:,7:end-0,:))/5;
% 
% 
% % 
% % corr image
% for roi=1:p.nROIs
% allCC{roi} = correlation_image_3D(b(:,:,:,roi), 4);
% end
% 
% % corr x mean image
for jj=1:p.nROIs
tmp = allCC{jj};
% tmp(tmp(:)>0)=sqrt(tmp(tmp(:)>0));
im_masks{jj} = Seg.mean_image{jj}.*tmp; 
end

%%

%Create masks
 im = im_masks;

    Seg.method = 'sobel';
    for roi=1:p.nROIs
    
    im{roi}(isnan(im{roi}))=0; 
    [~,threshold]=edge(im{roi}, Seg.method);    Seg.threshold(roi) = 0.5*threshold;
    BW      = edge(im{roi}, Seg.method, Seg.threshold(roi) );
    BWdil   = imdilate(BW, strel('square',3));
    Seg.fullmasks{roi} = imfill(BWdil, 'holes');%.*keeppix{roi};    Do this multiplication afterwards

    % keep maximal connected components
    Seg.masks{roi} = false(size(Seg.fullmasks{roi}));
    Seg.masks{roi}( get_bigg_conncomp(Seg.fullmasks{roi}) ) = true;
    Seg.masks{roi} = Seg.masks{roi} .*keeppix{roi};
    
    % Make background mask    
    Seg.all_background_masks{roi} =  ~(imdilate(Seg.fullmasks{roi}, strel('disk',4,4))); % inverse of detected ROI components
    Seg.all_background_masks{roi} = Seg.all_background_masks{roi} .*keeppix{roi};
    
    % Dark background
    bg_lim = prctile(Seg.mean_image{roi}(Seg.fullmasks{roi}),10);       % lowest 10th of fullmask
    Seg.background_masks{roi} = Seg.all_background_masks{roi} .*(Seg.mean_image{roi}<bg_lim);
    
    end   

    function pix = get_bigg_conncomp( bin_img )
components     = bwconncomp(bin_img);
numPixels      = cellfun(@numel,components.PixelIdxList);
[~,idx]        = max(numPixels);
pix            = components.PixelIdxList{idx(1)};
end