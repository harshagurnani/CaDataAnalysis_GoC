% Load DLC whisking angle and extract whisking variables
%
%
% Required input:
%    dlc_whisk_angle (whisker angle in radians)
%    fs (sampling rate in Hz )
% 
% Output:
%    dlc_whisk_time      time for whisking data
%    whisk_angle_filt    filtered whisking angle
%    whisk_set_point     whisker set point
%    whisk_amp           Whisking amplitude
%    whisk_phase  
% Modified by HG, from NACG's code.

function [whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase] = get_whisking_variables( dlc_whisk_angle, fs )

    

    % Get whisking resting point
    whisk_rest = median(dlc_whisk_angle);
    dlc_whisk_angle = dlc_whisk_angle - whisk_rest;

    % Filter to get denoised whisker angle
    fc = 30;
    [b,a] = butter(4,fc/(fs/2));
    whisk_angle_filt = filter(b,a,dlc_whisk_angle);

    whisk_angle_filt = filter(b,a,whisk_angle_filt(end:-1:1));

    whisk_angle_filt = whisk_angle_filt(end:-1:1);

    % To get set point, median filter by 500 ms
    %whisk_set_point = smoothdata(whisk_angle_filt,'movmedian',round(0.5* fs));
    whisk_set_point = smoothdata(whisk_angle_filt,'gaussian',round(0.5* fs));
    %smoothdata(whisk_angle_filt,'gaussian',[round(0.5* fs) 0] *2);

    % To get whisking amplitude and phase,
    fc2 = 30; fc1 = 8;
    [b,a] = butter(4,[fc1,fc2]/(fs/2));
    whisk_angle_bp = hilbert(filter(b,a,dlc_whisk_angle));
    
    whisk_phase = angle(whisk_angle_bp);
    whisk_amp = abs(whisk_angle_bp);
        
    % Convert to degrees
    whisk_set_point = rad2deg(whisk_set_point);
    whisk_amp = rad2deg(whisk_amp);

end