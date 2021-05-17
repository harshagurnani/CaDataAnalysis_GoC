function [lowerLid, upperLid, edist] = get_eyelid_dist( eye_vid, threshold, TN )
%% Compute Eyelid distance by computing median of each row to get a vector parallel to the minor axis of the eye. After thresholding, the transition from 1 to 0 to 1 gives poition of lower and upper eyelids

    edist = nan(TN,1);
    upperLid = edist;
    lowerLid = edist;
    
    initframe = imread( eye_vid, 1);
    if size(initframe,1) >size(initframe,2)
        initframe = transpose(initframe);
        do_tr = 1;
    else
        do_tr=0;
    end
    npixels = size(initframe, 2);
    width = size(initframe,1);
    clear initframe
    
    
    parfor frame = 1:TN
       im = imread(eye_vid, frame);
       minoraxis = nan(npixels,1);
       if size(im, 3) == 3
           im = rgb2gray(im);
       end
       if do_tr
           im =im';
       end
           
       for jj=1:npixels
           minoraxis(jj) = median( im(:,jj) );
       end
       minoraxis = minoraxis>threshold;
       ee = edge(minoraxis);
       ee = find(ee==1);
       if size(ee,1) == 1
           if ee > npixels/2
               ee = [1 ee];
           else ee = [ee npixels];
           end
       elseif size(ee,1) == 0
           ee = [1 npixels];
           
       end
       upperLid(frame) = ee(1);
       lowerLid(frame) = ee(end);
       
       edist(frame) = abs(lowerLid(frame) - upperLid(frame));
    end

    
end