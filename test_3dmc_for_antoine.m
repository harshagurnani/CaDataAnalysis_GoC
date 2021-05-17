nROI = numel(vol.start_pixels);

%Transformation matrix
AV = cell(nROI,1);

%Offset matrix
all_offsets = cell(nROI,1);

info = get_info_from_filename(parameters,'source',parameters.data_folder,'ROIs',parameters.ROIs,'repeats',parameters.repeats(1),'data_type',parameters.data_type);

%This section can be parallelized
for roi = nROI
    
v1 = [vol.start_pixels{roi}(:,1)-vol.stop_pixels{roi}(:,1)];
v2 = [vol.start_pixels{roi}(:,1)-vol.start_pixels{roi}(:,end)];
v3 = 0; 
v1 = v1./norm(v1);
v2 = v2./norm(v2);

% NOTE: no need for peakcorr for now (can check different weighting scheme later)
% [all_offsets{roi}, peakcorr{roi}] = get_offsets(info.files(roi,:), info.timepoints);
[all_offsets{roi}] = get_offsets(info.files(roi,:), info.timepoints);

% all offsets is the offsets in the plane - different axes for each patch - {roi}(2xt) -
all_offsets{roi} = permute(all_offsets{roi}, [2,1]);

% transform/ plane direction matrix
AV{roi} = [v1';v2'];    %should be 2x3

%
%    ---------   NOTE:  ---------------
% original code is solving an underdetermined system by assuming v3=0, but
% that is a WRONG assumption, as v3 is not necessarily aligned with z-axis!

% all_offsets{roi} = [ all_offsets{roi}(:,1) * v1(1) + all_offsets{roi}(:,2) * v2(1),...
% all_offsets{roi}(:,1) * v1(2) + all_offsets{roi}(:,2) * v2(2),...
% all_offsets{roi}(:,1) * v1(3) + all_offsets{roi}(:,2) * v2(3)];


end
% roi are concatenated:
all_offsets = cat(1, all_offsets{:} );  %reshaped to (2*roi)xt 
AV = cat(1, AV{:} );                    %reshaped to (2*roi)x3

nTimepoints = size(all_offsets,2);

%   ---------  NOTE:    ----------------
% System to solve is :
% At time(t): AV * o(t)' = all_offsets(:,t), where o is offset in cartesian
% coordinates. o(t) = (o_x(t), o_y(t), o_z(t) )

o = nan( nTimepoints, 3 );
% Can be parallelized - 
for t=1:nTimepoints
%   ---------NOTE:  ----------
% least squares fit - currently using default weight fun - need to think what works better
   o(t, : ) = tmp(2:4)';
end
all_offsets = o;    %now all_offsets is size t x 3 = (o_x,o_y,o_z) at each timepoint (row)


% ------ < OLD CODE > ------------
% all_offsets = cell2mat(reshape(all_offsets,[1,numel(all_offsets),1]));
% peakcorr = repmat(cell2mat( reshape(peakcorr, 1, 1, numel(vol.start_pixels)) ), [1, 3,1]);
% peakcorr = peakcorr./sum(peakcorr,3);
% all_offsets = nansum(all_offsets.*peakcorr,3);
% all_offsets = nanmean(all_offsets,3);
% ---------------------------------

% Calculate offsets in plane:
new_offsets = {};
for roi = 1:numel(vol.start_pixels)
v1 = [vol.start_pixels{roi}(:,1)-vol.stop_pixels{roi}(:,1)];
v2 = [vol.start_pixels{roi}(:,1)-vol.start_pixels{roi}(:,end)];
v1 = v1./norm(v1);
v2 = v2./norm(v2);
new_offsets{roi} = [ all_offsets * v1,...
                     all_offsets * v2];
end 


