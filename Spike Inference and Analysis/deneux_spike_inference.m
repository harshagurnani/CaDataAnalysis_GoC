function [ spk, fit, drift, parest ]  = deneux_spike_inference( x, dt, varargin )
%% Spike Inference using generative model from Deneux 2016
% Usage:
%   Inputs:
%       x:  Cell array OR (2D/3D)Matrix - Observed fluorescence
%       dt: Acquisition dt in ms
%   Output:
%       spk: Cell array of spike times for each roi
%       fit: Predicted fluor

%% Parse fluorescence input
if iscell(x)
    nROI = numel(x);
else
    % if x is a matrix
    if length( size(x) )==2
        nROI = size(x,2); x = mat2cell( x, size(x,1), ones(1, nROI) ); 
    elseif length( size(x) )==3
        nROI = size(x,3); x = mat2cell( x, size(x,1)*size(x,2), ones(1, nROI) );
    else
        error('Give x in valid format : cell array of 1xROI or a vector of [Time, (Trial), ROI]' )
    end
end

%% Set parameters for spike inference
[ par_all, redo_drift, driftrate_perc, drift_method, run_inference ] = parse_options( nROI, dt, varargin{:} );


% Set the range for drifting baseline
if redo_drift    
    for jj = 1:nROI,   par_all(jj).F0 = [ 0.5*prctile( x{jj}, 10) , 1.2*prctile( x{jj}, 80 ) ]; end
end

% Set the drift parameter
switch drift_method
    case 'fraction_dff'
        for jj = 1:nROI,   par_all(jj).drift.parameter = driftrate_perc * mean( par_all(jj). F0 ); end
        
    case 'identical'
        xmat = cell2mat(x);
        for jj = 1:nROI,   par_all(jj).drift.parameter = driftrate_perc * mean( xmat(:) ); end

end


%% Run algorithm from Deneux et al, 2016

switch run_inference
    case 'individual'
        % Calculate observation noise separately for each ROI
        disp( 'Calculating observation noise and running inference individually for each ROI...')
        for roi = 1:nROI,        roi
            [spk{1,roi}, fit{1,roi}, drift{1,roi}, parest{1,roi}]=spk_est( x(roi), par_all(roi));   end
    case 'batch'
        % same observation noise using power spectrum of all F
        disp( 'Calculating common observation noise for all ROI...')
        [spk, fit, drift, parest]=spk_est( x, par_all);
end


% for jj=1:nROI,  parest(jj).drift.method = drift_method;     end

end




%% parse optional arguments for inference parameters
function [ par_all, do_drift, driftrate, drift_method, run_inference ] = parse_options( nROI, dt, varargin )
% Set parameters for spike inference

    par = spk_est('par'); par.dt = dt;
    par_all(1:nROI) = deal( par );
    % default values
    do_drift = true; 
    driftrate = 0.01;
    drift_method ='fraction_dff';
    estimate = 'map';
    dorise = false;
    rise_time = 15; %ms
    run_inference = 'batch';
    
    dorise_reset = false;
    nargs = length(varargin);

    if nargs>1
        % check for par structure
        
        for jj=1:nargs-1
            if strcmp(varargin{jj}, 'par'),     par = varargin{jj+1};       par.dt = dt;              par_all(1:nROI) = deal( par );end
            if strcmp(varargin{jj}, 'par_all'), par_all = varargin{jj+1};   par_all(:).dt = deal(dt); end
        end
        jj=1; 
        while jj<=nargs
            
            if strcmp( varargin{jj}, 'dodrift')
                do_drift = varargin{jj+1};
                jj=jj+2;
            elseif strcmp( varargin{jj}, 'driftrate')
                driftrate = varargin{jj+1};
                jj=jj+2;
            elseif strcmp( varargin{jj}, 'drift_method')
                drift_method = varargin{jj+1};
                jj=jj+2;    
            elseif strcmp( varargin{jj}, 'est_method')
                estimate = varargin{jj+1};
                jj=jj+2;
            elseif strcmp( varargin{jj}, 'dorise')
                dorise = varargin{jj+1};
                jj=jj+2;
            elseif strcmp( varargin{jj}, 'rise_time')
                dorise_reset = true;
                rise_time = varargin{jj+1};
                jj=jj+2;
            elseif strcmp( varargin{jj}, 'run_method')
                run_inference = varargin{jj+1};
                jj=jj+2;
            else
                error( [varargin{jj}, ':  Invalid argument name'] )
            end
        end
    end
    
    if dorise_reset,        dorise = true; end
    

    
    % Set the estimation algorithm
    %Use MAP for a single spike train - the map estimate. or 'Samples' for returning multiple spike trains?
    for jj =1:nROI,  par_all(jj).algo.estimate = estimate;  end

    if dorise
        % Use rise time for GCamp6f
        for jj = 1:nROI,   par_all(jj).ton = rise_time/1000; end
    end
    
    if (~do_drift && isscalar(par_all(1).F0) ),   driftrate = 0; end

end

