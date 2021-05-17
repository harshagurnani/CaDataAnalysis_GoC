% speed = importdata('Speed Data 001.txt');
absTime = speed.textdata(:,2); 
relTime = speed.data(:,1)*.001; %in ms

% convert absTime from string to num
absTime = arrayfun( @(jj) ( str2double(absTime{jj}(1:2))*60*60 + str2double(absTime{jj}(4:5))*60 + str2double(absTime{jj}(7:end)) )*1000, 1:length(absTime));% in ms
absTime = absTime';

% find encoder resets - start and end of trial
flips = find(relTime(2:end)<relTime(1:end-1))+1;    

exp_start = 1;  % experimental start trigger
% check for delayed trigger?
tlen = prctile(diff(flips), 80 ); 
if relTime(flips(1)-1) < 0.5*tlen, exp_start = flips(1); flips = flips(2:end);   end
absTime = absTime - absTime(exp_start) + relTime(flips(1)); %set exp_start trigger time as absolute zero 
%     
%%%% OR comment out to load Speed abs time, reltimes, MI
%%
% Trial start and end time - absolute times
nTrials = floor(length(flips)/2)+1;
TrialTime = nan( nTrials, 2 );
TrialTime(1,1) = 0; TrialTime(1,2) = absTime( flips(1) ); 
for jj=2:nTrials
    TrialTime(jj,1) = absTime(flips(2*(jj-1)))-relTime(flips(2*(jj-1))); % Trial start time = abstime at flip - rel time at flip
    if jj==nTrials, TrialTime(jj,2) = absTime(end); 
    else, TrialTime(jj,2) = absTime(flips(2*jj-1));                      % trial end time = abstime at next flip
    end
end

%%
% Load MI into workspace // Preprocess your MI 
load('video.mat','MI');    % First column is MI, second column is time

TrialMI = cell( nTrials, 1); 
for jj=1:nTrials
    si = find( MI(:,2) > TrialTime(jj,1), 1);   si = max(1, si-1) ;
    li = find( MI(:,2) > TrialTime(jj,2), 1);   if isempty(li), li = length(MI); end
    TrialMI{jj} = MI( si:li, :);
    TrialMI{jj}(:,2) = TrialMI{jj}(:, 2)-TrialTime(jj,1);   % subtract absolute trial start time to convert to relative time within trial
end

%% Plotting
figure; hold on; cm=colormap(redblue(nTrials));
for jj=1:nTrials, plot( TrialMI{jj}(:,2), TrialMI{jj}(:,1), 'color', cm(jj,:)); end
xlim([0 TrialMI{jj}(end,2)])