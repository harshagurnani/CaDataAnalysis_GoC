function [FEC] = get_fec( eye_vid,  filt_size, threshold, TN )
%% Compute Fractional Eyelid Closure (FEC) by median filtering, thresholding (to binary) and getting a normalized/fractional sum such that FULL CLOSURE FEC = 1 and FULLY OPENED FEC = 0

    FEC = nan(TN,1);
    if size(filt_size, 2) == 1
        filt_size = [filt_size filt_size];
    end
    initframe = imread( eye_vid, 1);
    npixels = size(initframe, 1) * size(initframe,2) ;
    clear initframe
    parfor frame = 1:TN
       im = imread(eye_vid, frame);
       if size(im, 3) == 3
           im = rgb2gray(im);
       end
       im = medfilt2( im, filt_size );
       FEC(frame) = sum(im(:) > threshold);
    end
       
    m_fec = min(FEC);
    FEC = (FEC-m_fec)./(npixels - m_fec);
    
end