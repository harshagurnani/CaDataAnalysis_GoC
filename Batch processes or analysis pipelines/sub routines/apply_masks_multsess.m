[ Seg.(sess).r,  Seg.(sess).t ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.masks, 'roi_indx', Seg.mask_roi);         % avg somatic fluorescence
[ Seg.(sess).background, ~ ] = mask_nonbinary_and_avg( Basic.(sess).r, Basic.(sess).t, Seg.background_masks);                           % avg dark baclground fluorescence

% remove uncorrected frames
errorid = find(MC.(sess).error);
Seg.(sess).error = MC.(sess).error;
if max(diff(errorid))<8, Seg.(sess).error(errorid(1):errorid(end)) = true; end
Seg.(sess).r( Seg.(sess).error,: ) = NaN;


