%% P5 - Compute correlations

% PAIRWISE

Corr.pw.params.ds_rate             = 20;           % Hz
Corr.pw.params.maxlag              = 2000;         % ms
Corr.pw.params.nshuffle            = 500;
Corr.pw.params.loco.extratime      = 500;          % ms
Corr.pw.params.whisk.extratime     = 500;          % ms
Corr.pw.params.puff.extratime      = [200 1000];   % ms
Corr.pw.params.loco_maxlag         = 2000;         % ms
Corr.pw.params.whisk_maxlag        = 2000;         % ms
Corr.pw.params.puff_maxlag         = 2000;         % ms
Corr.pw.params.whisk_doconcat       = false;
Corr.pw.params.loco_doconcat       = false;
Corr.pw.params.interSess           = 2000;         % ms

% WITH BEHAVIOUR
Corr.beh.params.ds_rate             = 10;           % Hz
Corr.beh.params.maxlag              = 2000;         % ms
Corr.beh.params.nshuffle            = 500;


%% %---------------------------------------------------
if get_corr
    get_correlations_for_dff;
end 

save([exp.analysed, '\correlations.mat'], 'Corr','p','-v7.3');
disp('.....................Saved Pairwise and Behavioural Correlations')
%--------------------------------------------------



%% P7 - Determine Airpuff Responses

if get_puff
    PuffResponses.params.Window = [-550 1550]; %ms                             % Start and end Time around airpuff event to crop (0 = airpuff)
    PuffResponses.params.baselineWindow = [-500 -300]; %ms                     % for subtracting mean baseline (pre-event) activity

    PuffResponses.params.dsRate = 20;  %Hz                                     % Sampling rate for averaged responses aligned to puff
    PuffResponses.params.smoothLoco  = 100; %ms
    PuffResponses.params.smoothWhisk = 100; %ms

    PuffResponses.params.noLoco_MeanE = 0.1;
    airpuff_analysis 
    
end

%% P8 - Dimensionality reduction
PC_Corr.params.ds_rate  = 20;
PC_Corr.params.maxlag   = 2000;
PC_Corr.params.nshuffle = 500;

run_pca_dff;
