function [ best_traj, R, MC_correct ] = complete_mc( RawData, maxD, lambda, varargin )

tic

%Max time and number of Patches
Tm = size(RawData,3);
nTrials = size(RawData,4);
nROI = size(RawData,5);
nX = size(RawData,2);
nY = size(RawData,1);

% DX and DY are the maximum x and y displacements. Row number is x, column
% number is y. 
if size(maxD,2) == 1
    DX = maxD;
    DY = maxD;
elseif size(maxD,2) == 2
    DX = maxD(2);
    DY = maxD(1);
end
del_x = -DX:DX;
del_y = -DY:DY;

%Frame to analyse: Fy x Fx
Fx = 1+DX:nX-DX;
Fy = 1+DY:nY-DY;

%% Parsing optional arguments
gamma=[];
montage = [];
change_ref = [];
refscale=[];
if nargin > 3  

    charctr=0;
    numctr=0;
    start_read_var = 4;
    while start_read_var <= nargin
        vn = start_read_var-3;
       if ischar( varargin{vn} )
           charctr=1;
           if strcmp(varargin{vn}, 'gamma' )
               gamma = varargin{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargin{vn}, 'montage' )
               montage = varargin{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargin{vn}, 'change_ref' )
               changeref = varargin{vn+1};
               start_read_var = start_read_var + 2; 
           elseif strcmp(varargin{vn}, 'refscale' )
               refscale = varargin{vn+1};
               start_read_var = start_read_var + 2;    
           else 
                error('Argument name %s is not an optional argument', varargin{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   gamma = varargin{vn};
                   start_read_var = start_read_var+1;
               case 2
                   montage = varargin{vn};
                   start_read_var = start_read_var+1;
               case 3
                   changeref = varargin{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end

if isempty(gamma)
    gamma = 10;
end

if isempty(montage)
    montage=false;
end

if isempty(changeref)
    changeref=false;
end

if isempty(refscale)
    refscale=30;
end
%------------------------------------------------------------------------%
%       Transition matrix for displacements delta_x, delta_y
%------------------------------------------------------------------------%
% T( (dx, dy) --> (dx',dy') ) = 1/(2*pi*lamda^2) * exp(-r/lambda) where
% lambda is spread of probablitity fn and r is Euclidean distance between
% successive displacements i.e. r = sqrt( (dx'-dx)^2 + (dy'-dy)^2 )

x2 = repmat(del_x, 2*DX + 1, 1) - repmat(del_x', 1, 2*DX +1);
x2=x2.^2;
y2 = repmat(del_y, 2*DY + 1, 1) - repmat(del_y', 1, 2*DY +1);
y2= y2.^2;
r = sqrt( repmat( y2, [1,1, 2*DX+1, 2*DX+1] ) + ...
          repmat(reshape(x2, [1,1,2*DX+1, 2*DX+1]), [2*DY+1,2*DY+1,1,1]) );

% No need to calculate probabilities as we are using log Prob
% T_del = 1/(2*pi*lambda^2) * exp( - r/lambda );
% T_del = permute(T_del, [1,3,2,4]);

% for kk=1:2*DX+1
%     for jj=1:2*DY+1
%         T_del(jj,kk,:,:) = T_del(jj,kk,:,:)./sum(sum(T_del(jj,kk,:,:)));
%     end
% end


%--------------------------
% Photon count parameter: gamma = 1 if ref R and obs idata I are photon
% counts, otherwise use calibration factor from intensity to counts


%%
% Reference Image
R = nan(nY, nX, nROI);

% Reference is plane with minimum squared distance from next frame? 
% RawData : Y x X x Time x ROI
parfor roi = 1:nROI
        mi_roi = simple_motion_index( squeeze(reshape(RawData(:,:, :, :, roi), [nY, nX, Tm*nTrials] )), repmat([1:Tm*nTrials]', [1,2]) );
        M = max(mi_roi(:,1));
        indx = find( mi_roi(:,1) == M );
        t = mod(indx, Tm)+1;
        trial = ceil(indx/Tm);
        R(:,:, roi) = squeeze(RawData(:,:, t, trial, roi));
end
R = R(Fy,Fx,:);


% 
% lambda = [0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2];
%  gamma = [1, 5, 10, 20, 50, 100];

% lambda = [0.1,  2];

% maxLP = -Inf;
% best_lm = [];


maxlogP = nan( nTrials, 1 );
MC_correct = nan( size(Fy,2), size(Fx,2), Tm, nTrials, nROI);
if montage
    best_traj = nan( nTrials, Tm ); 
else
    best_traj = nan( nTrials, Tm*nROI );
end
parfor trial=1:nTrials
    
%     [ best_traj(trial,:), maxlogP(trial) ]  = do_mc_use_viterbi( DX, DY, r, 'uniform', squeeze(floor(RawData(:,:,:,trial,:)*10)/10), squeeze(floor(R(:,:,trial,:)*10)/10),  Tm, nROI, lambda, gamma );
    if montage
       if changeref
        [ best_traj(trial,:), maxlogP(trial) ]  = do_mc_change_ref_montage( DX, DY, r, 'uniform', squeeze(RawData(:,:,:,trial,:)), R, refscale, Tm, nROI, lambda, gamma );
       else
        [ best_traj(trial,:), maxlogP(trial) ]  = do_mc_use_viterbi_montage( DX, DY, r, 'uniform', squeeze(RawData(:,:,:,trial,:)), R,  Tm, nROI, lambda, gamma );
       end
    else
       if changeref
        [ best_traj(trial,:), maxlogP(trial) ]  = do_mc_change_ref2( DX, DY, r, 'uniform', squeeze(RawData(:,:,:,trial,:)), R, refscale,  Tm, nROI, lambda, gamma );
       else
        [ best_traj(trial,:), maxlogP(trial) ]  = do_mc_use_viterbi( DX, DY, r, 'uniform', squeeze(RawData(:,:,:,trial,:)), R,  Tm, nROI, lambda, gamma );
       end
    end
end
% disp_x =  ceil( best_traj(:,(0:Tm-1)*nROI + roi)/(2*DY+1) ) - (DY+1);
% disp_y =  mod( best_traj(:,(0:Tm-1)*nROI + roi) -1, 2*DY+1 ) - DY;
% 
% parfor roi = 1:nROI
%     temp = RawData(:,:,:,:,roi);
%     for t=1:Tm
%     mc = (t- 1)*nROI + roi;
%     for trial = 1:nTrials
%         MC_correct(:,:,t,trial, roi) = temp(Fy+disp_y(trial, mc), Fx+disp_x(trial, mc), t, trial );
%     end
%     end
% end
% for lm = 1:size(lambda,2)    
%     if maxlogP(lm) > maxLP
%         maxLP = maxlogP(lm);
%         best_lm = lambda(lm);
%     end
% end
disp_x =  ceil( best_traj/(2*DY+1) ) - (DX+1);
disp_y =  mod( best_traj -1, 2*DY+1 ) - DY;
nfx = size(Fx,2);
nfy = size(Fy,2);
if montage

    for trial = 1:nTrials
        for t=1:Tm
            for roi=1:nROI
             MC_correct(:,:,t, trial, roi) = RawData( Fy + disp_y(trial,t), Fx + disp_x(trial,t), t, trial , roi);
            end
        end
    end
else
for roi = 1:nROI
    disp_x =  ceil( best_traj(:,(0:Tm-1)*nROI + roi)/(2*DY+1) ) - (DX+1);
    disp_y =  mod( best_traj(:,(0:Tm-1)*nROI + roi) -1, 2*DY+1 ) - DY;
    temp = RawData(:,:,:,:,roi);
    mc_roi = nan(nfy, nfx, Tm, nTrials);
    for t=1:Tm
        for trial = 1:nTrials
            mc_roi(:,:,t, trial) = temp( Fy + disp_y(trial,t), Fx + disp_x(trial,t), t, trial );
        end
    end
    MC_correct(:,:,:,:, roi) = mc_roi;
end
end

totaltime = toc
end
% 
% 
% 
% %% Viterbi Algorithm
% 
% function [ best_traj, maxlogP ] = do_mc_use_viterbi( DX, DY, r, initial_st, I, R, Tm, nROI, lambda, gamma )
% 
% %-------------------------------------------------------------------------%
% %                                PARAMETERS
% %-------------------------------------------------------------------------%
% 
% % Patch size
% nX = size(I,2);
% nY = size(I,1);
% 
% % DX and DY are the maximum x and y displacements. Row number is x, column
% % number is y. 
% del_x = -DX:DX;
% del_y = -DY:DY;
% 
% %Frame to analyse: Fy x Fx
% Fx = 1+DX:nX-DX;
% Fy = 1+DY:nY-DY;
% 
% %-------------------------------------------------------------------------%
% %                          INITIALIZATION
% %-------------------------------------------------------------------------%
% 
% % Initialise to uniform prob [or set probability = 1 for state that
% % maximizes 2-D corr for first ROI.]
% if strcmp( initial_st, 'uniform' ) 
%     initial_state_prob = ones(2*DY+1,2*DX+1) * 1/(2*DY+1) * 1/(2*DX+1);   
% elseif strcmp( initial_st, 'max-2D')
%     % Use xcorr2 to align first ROI with its reference     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% 
% prev_state_prob = initial_state_prob;
% curr_state_prob = prev_state_prob;
% 
% % Trajectory until t-1 or t. Preallocate for speed,
% % Values are linear indices corresponding to displacement pair at each
% % MC_cycle.
% traj_prev = nan( (2*DY+1)* (2*DX+1), Tm*nROI );
% traj_curr = traj_prev;
% 
% % Probability of the trajectories until elast completed MC_cycle. To find
% % best trajectory in the end.
% % P( Traj_i,j,t ) = P( del_t = i,j | Traj_i',j',t-1) * P(Traj_i',t',t-1) 
% % log traj_prob( i,j,t) = log P (i,j | i',j', Traj_i', j', t-1) + log P(i',j',t-1)
% traj_prob = zeros( (2*DY+1)* (2*DX+1),1);
% traj_prob_new = traj_prob;
% 
% % Linear indexing of all possible displacement pairs
% del_index = reshape(1:(2*DX+1)*(2*DY+1), [(2*DY+1), (2*DX+1)]);
% 
% %------------------------------------------------------------------------%
% %               MOTION CORRECTION: Each Viterbi Probability
% %------------------------------------------------------------------------%
% tic
% 
% for MC_cycle = 1:Tm*nROI
%     
%    t = ceil( MC_cycle/nROI);
%    plane = mod(MC_cycle-1, nROI) + 1;
%    
%    % Fo each possible displacement( dy(t), dx(t) ):
%    for dx_t = 1:2*DX+1
%    for dy_t = 1:2*DY+1
%        
%         maxlogP = -Inf;
%         del_tminus1 = nan;
%         
%         logP1 = sum(sum( probIgivenR( I(Fy+ del_y(dy_t), Fx+ del_x(dx_t), t, plane), R(:,:, plane ), gamma) ) );
%         
%         for dx_tminus1 = 1:2*DX+1
%         for dy_tminus1 = 1:2*DY+1
%             
%        % Compute (log) Probability for current displacement = (i, j) for
%        % each displacement pair (i', j') in previous step. Choose
%        % displacement at t-1 [del(t-1)] for del(t) = (i,j) corresponding to
%        % the evaluated current state that maximizes probability of
%        % trajectory passing through (i,j) at time t. This is the 
%        % log viterbi_traj_prob( i,j,t) 	i.e. v(i,j,t)
%        
%        % v(i,j,t) = log P(i, j, Traj_i,j,t)
%        % v(i,j,t) = max_{i',j'}     < log P (i,j, I | i',j', Traj_i',j',t-1, R ) + log P( i', j', Traj_i',j',t-1) >
%        %          = max_{i',j'}     < log P(i,j, I | i',j', R)   [MARKOV]        + v(i',j',t-1) >
%        %          = max_{i',j'}     < log P( I(t)| R, i,j ) + log T(i',j'-->i,j) + v(i',j',t-1)  >
%        
%        % Term 1:    log P(I(t) | R, i, j) = sum_{x,y} [ log P ( ObsI(x,y,plane) | RefR(x+i, y+j, plane) ) ]
%        % Term 2:    log T(i', j' --> i,j)
%        % Term 3:    v(i',j',t-1)
%        
%            logP2 = -r( dy_tminus1, dx_tminus1, dy_t, dx_t )/lambda; 
%            logP = logP1 + logP2 + traj_prob( del_index(dy_tminus1, dx_tminus1) );
%            
%            if logP > maxlogP
%                maxlogP = logP;
%                del_tminus1 = del_index(dy_tminus1, dx_tminus1);
%            elseif logP == maxlogP
%                del_tminus1 = [del_tminus1, del_index(dy_tminus1, dx_tminus1)];
%            end
%            
%         end
%         end
%         
%         %%% If multiple del(t-1) give same logP
%         if size(del_tminus1,2) > 1
%            % choose which one maximizes 2D-corr      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            del_tminus1 = del_tminus1(1);
%         end
%         
%         traj_prob_new( del_index(dy_t, dx_t) ) = maxlogP;
%         traj_curr( del_index(dy_t, dx_t), 1:MC_cycle ) = cat(2, traj_prev(del_tminus1, 1:MC_cycle-1), del_index(dy_t, dx_t) ); 
%         
%    end
%    end
%    
%    traj_prob = traj_prob_new;           % Update and set to v(i,j,t)
%    traj_prev(:,1:MC_cycle) = traj_curr(:,1:MC_cycle); % Update to trajectory upto t-1
%     
% end
% 
% %------------------------------------------------------------------------%
% %           Find Best Trajectory
% %------------------------------------------------------------------------%
% maxlogP = max(traj_prob);
% best_traj = traj_curr( (traj_prob == maxlogP), 1:Tm*nROI );
% 
% maxlogP = maxlogP - Tm*nROI*log(2*pi*lambda*lambda);  %Add the lambda-dependent contribution
% 
% toc
% 
% end
% 
% %% Probability Function
% function logprobs = probIgivenR( ObsI,  MCRef, gamma)
% 
%     ObsI = squeeze(ObsI);   % Observed intensity
%     MCRef = squeeze(MCRef); % Correct sub-frame after MC - "Real mean for ObsI"
%     
%     % For each pixel in plane:
%     % P(I_x,y,p | R_x,y,p) = (gamma*R)^(gamma*I) exp(-gamma*R)  / (gamma*I)!
%     % log P(I|R) = gamma*( I * log R - R) + constant K(gamma, I)
%     logprobs = gamma * ( ObsI.* log(MCRef) - MCRef );
% 
% end