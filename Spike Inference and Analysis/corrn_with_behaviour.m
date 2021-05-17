

%% Neural responses
% Raw data present?
PCOR.do_raw = false;
if exist('data','var'), PCOR.do_raw = true; end

% Load spike trains and raw trains
disp('Loading data and spike trains')
if exist('spk','var')
    FR = spk;       FR_t = xt;
    if iscell(FR{1}),for jj=1:numel(FR), FR{jj}=FR{jj}{1}; end; end
    if exist('r','var'), raw = r; raw_t = t; PCOR.do_raw = true; end
elseif exist('spk_pooled', 'var')
   FR = spk_pooled; FR_t = pooled_t;
   if iscell(FR{1}),for jj=1:numel(FR), FR{jj}=FR{jj}{1}; end; end
   if PCOR.do_raw, raw = data.group.raw; raw_t = data.group.time; end
elseif exist('spk_ind', 'var')
    FR = spk_ind;   FR_t = ind_t;
    if PCOR.do_raw, raw = data.ind.raw; raw_t = data.ind.time; end
end

if iscell(FR),   FR   = cell2mat(FR);   end
if iscell(FR_t), FR_t = cell2mat(FR_t); end
if PCOR.do_raw, if iscell(raw), raw = cell2mat(raw); end ; end
FR = double(FR);
POI = [0 max(FR_t(:))];    %over entire timeseries

%% Normalise raw data to get dFF
disp('Normalising raw data...')
if PCOR.do_raw
    if length(size(raw))<3
        warning('Reshaping your raw data!')
        raw   = reshape( raw,   [p.nTimepoints, length(p.trials), size(raw,2)] );
        raw_t = reshape( raw_t, [p.nTimepoints, length(p.trials), size(raw_t,2)] );
        
    end
    [norm_data, baseline] = normalize_and_smooth(raw, raw_t, p); 
    norm_data = reshape(norm_data, [p.nTimepoints* length(p.trials), size(norm_data,3)]);
end
    
%Smooth spike train into instantaneous firing rates
PCOR.smooth_bin.nrn = 100;%ms
dt = 1000/p.acquisition_rate; sfactor = 2*floor(0.5*PCOR.smooth_bin.nrn/dt)+1;
for jj=1:size(FR,2), FR(:,jj) = 1000/dt*smooth(FR(:,jj), sfactor );    end

%% PC projections
%check if there is a projection into PCA subspace and smooth it
do_PC = false;
if exist( 'PCA_result', 'var')
    do_PC = true;   disp('Loaded PCA results')
    PC_comp = PCA_result.proj_orig'; ndimPC = size(PC_comp,2);  
%     for jj=1:ndimPC, PC_comp(:,jj) = 1000/dt*smooth(PC_comp(:,jj), sfactor);    end
else
    try
        disp('Running PCA on population firing rates')
        PCA_result = test_pca_on_fr( FR, FR_t);
        PC_comp = PCA_result.proj_orig'; ndimPC = size(PC_comp,2);  
%         for jj=1:ndimPC, PC_comp(:,jj) = 1000/dt*smooth(PC_comp(:,jj), sfactor);    end
        do_PC=true;
    catch
        disp('Could not run PCA on firing rates.. Proceeding without it')
    end
end

%% Behavioural variables
if     exist('speed', 'var'), speed_new = speed;
elseif exist('spd',   'var'), speed_new = spd;      end

if iscell(speed_new), speed_new = cell2mat(speed_new); end
speed_new(:,2)=abs(speed_new(:,2));
%Smooth speed 
PCOR.smooth_bin.speed = 40;%ms
spd_dt     = 2; %ms
sfactor = 2*floor(0.5*PCOR.smooth_bin.speed/spd_dt)+1;
speed_new(:,2) = sfactor*60/50*smooth( speed_new(:,2), sfactor );   %cm/s



%% Calculate correlations
% Correlation with running speed
disp('Calculating correlations with running speed... ZZZ ...')
PCOR.dsrate.speed = 10;   %Hz
PCOR.nshuffle.speed    = 500;
[ allcorr.fr.run.corr, allcorr.fr.run.sig, allcorr.fr.run.rho, allcorr.fr.run.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, speed_new(:,2), speed_new(:,1), ...
                                                                                                         FR, FR_t, PCOR.dsrate.speed, PCOR.nshuffle.speed, 0 );
if do_PC
[ allcorr.PC.run.corr, allcorr.PC.run.sig, allcorr.PC.run.rho, allcorr.PC.run.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, speed_new(:,2), speed_new(:,1), ...
                                                                                                     PC_comp, FR_t, PCOR.dsrate.speed, PCOR.nshuffle.speed, 0 );
end 

if PCOR.do_raw
[ allcorr.dff.run.corr, allcorr.dff.run.sig, allcorr.dff.run.rho, allcorr.dff.run.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, speed_new(:,2), speed_new(:,1), ...
                                                                                                     norm_data, FR_t, PCOR.dsrate.speed, PCOR.nshuffle.speed, 0 );
end 

% Correlation with motion index
disp('Calculating correlations with whisker pad motion index')
PCOR.dsrate.whisk = 10;   %Hz
PCOR.nshuffle.whisk    = 500;                                                                                                     
[ allcorr.fr.whisk.corr, allcorr.fr.whisk.sig, allcorr.fr.whisk.rho, allcorr.fr.whisk.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, mi_eye(:,1), mi_eye(:,2), ...
                                                                                                         FR, FR_t, PCOR.dsrate.whisk, PCOR.nshuffle.whisk, 0 );
if do_PC
[ allcorr.PC.whisk.corr, allcorr.PC.whisk.sig, allcorr.PC.whisk.rho, allcorr.PC.whisk.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, mi_eye(:,1), mi_eye(:,2),PC_comp, FR_t, PCOR.dsrate.whisk, PCOR.nshuffle.whisk, 0 );
end

if do_PC
[ allcorr.dff.whisk.corr, allcorr.dff.whisk.sig, allcorr.dff.whisk.rho, allcorr.dff.whisk.shuffcorr, ~] =  all_period_shuffle_zcorr( POI, mi_eye(:,1), mi_eye(:,2), ...
                                                                                                     norm_data, FR_t, PCOR.dsrate.whisk, PCOR.nshuffle.whisk, 0 );
end