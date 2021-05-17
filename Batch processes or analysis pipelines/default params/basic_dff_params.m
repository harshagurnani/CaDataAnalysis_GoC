function newBasic = basic_dff_params( oldBasic )
Basic.ROI_type            = 'patch'; % or 'points', 'plane'                 % What kind of data did you record?

Basic.ROIs                = 0;                                              % which ROI to keep for analysis;    0 = All, else give array with ROI numbers
Basic.Groups              = 0;                                              % Which groups to keep for analysis; 0 = all
Basic.trials              = 0;                                              % Which trials to keep for analysis; 0 = all
Basic.continuous_time     = true;                                           % Single continuous time - computed based on cycle lengths* number of timepoints
                                                                            % inter-trial time can be added if speed file is available.
                                                                            % If false, each trial is considered a repetition
% motion correction (for patches and plane timeseries)
Basic.MC                  = true;                                           % If true, corr-based motion correction will be applied to planes/patches before any segmentation/normalisation
Basic.MC_select_ROI       = true;                                           % If true, user input is required at the time of execution to select ROI used for MC. For patches, that is patch#, for planes, user can draw rectangles.
                                                                            % Can be set to an array (for patch#) or 5-column matrix (each row is a rectangle - [plane #, x1, xn, y1, yn]) instead of boolean value
                                                                            % If false, all patches are used for MC - can be memory intensive!!
Basic.MC_max_disp         = 8;                                              % Max displacement allowed in any MC-frame; If optimal disp> max_disp, F = NaN in that frame                                                                            


% segmentation/averaging (for patches and plane timeseries)
Basic.patch_avg           = false;                                           % If true, mask(s) is applied to patch, else each pixel is normalised independently.
Basic.patch_multiple      = false;                                          % If true, multiple masks are kept for a single patch(plane)

Basic.get_background      = true;                                           % Identify background_per_patch.
Basic.ROI_background      = [];                                             % Array of ROI numbers that were put as backgrounds

% normalisation and smoothing
Basic.norm_baseline       = 'percentile';    % or 'interval'                % How to compute F0?
Basic.base_interval       = [];%ms                                          % Time interval considered as baseline, and its time-averaged F = F0
Basic.base_percentile     = 10;                                             % percentile to fet F0 per trace
Basic.base_invert         = false;                                          % If true, events are assumed to be negative, and hence normalisation is changed accordingly. 
                                                                            % Top base_percentile rather than lowest base_percentile is considered as baseline
Basic.baseline_pertrial   = false;                                          % if true, each trial has its own baseline

Basic.smooth_scale        = 100;%ms                                         % Filter size for smoothing normalised traces (boxcar_filter)
Basic.smooth_causal       = false;                                          % If true, filter is causal i.e [t-n : t] rather than [t-n/2 : t+n/2]

% merge defaults into struct
allfields = fieldnames( Basic );
for ff = 1:numel(allfields )
   if ~isfield( oldBasic, allfields{ff} )
       oldBasic.(allfields{ff}) = Basic.(allfields{ff}); 
   end
end

newBasic = oldBasic;
end