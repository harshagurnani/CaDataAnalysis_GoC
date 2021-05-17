% Matlab Read Trial Time from Encoder

nTrials = length(s);
TrialTime = nan( nTrials, 2 );
for jj=1:nTrials
    TrialTime(jj,1) = s{jj}( 1,1)*1000 - s{jj}(1,2); %ms
    TrialTime(jj,2) = s{jj}(end,1)*1000; %ms
end


% TrialTime(1,1) = 0; TrialTime(1,2) = absTime( flips(1) ); 
% for jj=2:nTrials
%     TrialTime(jj,1) = absTime(flips(2*(jj-1)))-relTime(flips(2*(jj-1))); % Trial start time = abstime at flip - rel time at flip
%     if jj==nTrials, TrialTime(jj,2) = absTime(end); 
%     else, TrialTime(jj,2) = absTime(flips(2*jj-1));                      % trial end time = abstime at next flip
%     end
% end