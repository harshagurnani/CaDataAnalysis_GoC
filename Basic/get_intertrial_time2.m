function [trialtime, itt, varargout ] = get_intertrial_time( params, varargin )
%% get fifo time of all Trial trigger (Start and End of trial)length of trials and intertrial duration based on fifo trigger 
% for older data, get trial times and intertrial dur from encoder data

    exppath = params.exp_path;
    %---------------- compatibility with old data
    if ~isfield( params, 'maxTrial') 
        try params.maxTrial = params.maxTrials; catch, params.maxTrial = params.nTrials; end
    end
    allTrigger = [];
    
    %rewrite based on software version - check with sameer - qqqqq
    inifile = [exppath,'/Experiment Header.ini'];
    fid = fopen(inifile,'r');
    text = textscan(fid, '%s','Delimiter','\t');
    fclose(fid);
    text = text{1};

    %find section line of FIFO times
    fifoLine = find(contains(text, '[Intertrial FIFO Times]'));
    
    if ~isempty( fifoLine)
        try 
            allLines = text(fifoLine+1:fifoLine+(2*2*params.maxTrial) );
        catch
            allLines = [];
        end
    else
        allLines = [];
    end
    if ~isempty(allLines)
        %new version
        allTrigger = arrayfun(@(jj) str2double(allLines{jj}), 1:length(allLines));
        allTrigger = (reshape(allTrigger,[2, 2*params.maxTrial]))';
        allTrigger(:,2) = allTrigger(:,2)*1e3; %ms
        
        trialtime = reshape(allTrigger(:,2),[2,params.maxTrial])';
        itt = arrayfun( @(t) trialtime(t+1,1)-trialtime(t,2), 1:params.maxTrial-1 )'; 
        trialtime = trialtime( params.trials, :);   %Crop to requested trials
        
    else
       %read intertrial from encoder (encoder_shift_dt is added to
       %intertrial_duration: itt)
       [trialtime, itt] = trialtime_from_encoder2( params );

    end
    
    if nargout>2
        varargout{1} = allTrigger;
    end
    
end