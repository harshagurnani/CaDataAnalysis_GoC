function logprobs = probIgivenR( ObsI,  MCRef, gamma)
%% Probability Function

    ObsI = squeeze(ObsI);   % Observed intensity at MCcorrect-shifted subframe
    MCRef = squeeze(MCRef); % reference intensity at centered subframe
    
    % For each pixel in plane:
    % P(I_x,y,p | R_x,y,p) = (gamma*R)^(gamma*I) exp(-gamma*R)  / (gamma*I)!
    % log P(I|R) = gamma*( I * log R - R) + constant K(gamma, I)
    logprobs = gamma * ( ObsI.* log(MCRef) - MCRef );

end