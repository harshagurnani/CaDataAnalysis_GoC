function [ best_traj, maxlogP ] = do_mc_change_ref2( DX, DY, r, initial_st, I, R, refscale, Tm, nROI, lambda, gamma )
%% Viterbi Algorithm

%-------------------------------------------------------------------------%
%                                PARAMETERS
%-------------------------------------------------------------------------%

% Patch size
nX = size(I,2);
nY = size(I,1);

% DX and DY are the maximum x and y displacements. Row number is x, column
% number is y. 
del_x = -DX:DX;
del_y = -DY:DY;

%Frame to analyse: Fy x Fx
Fx = 1+DX:nX-DX;
Fy = 1+DY:nY-DY;

%-------------------------------------------------------------------------%
%                          INITIALIZATION
%-------------------------------------------------------------------------%

% Initialise to uniform prob [or set probability = 1 for state that
% maximizes 2-D corr for first ROI.]
if strcmp( initial_st, 'uniform' ) 
    initial_state_prob = ones(2*DY+1,2*DX+1) * 1/(2*DY+1) * 1/(2*DX+1);   
elseif strcmp( initial_st, 'max-2D')
    % Use xcorr2 to align first ROI with its reference     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

prev_state_prob = initial_state_prob;
curr_state_prob = prev_state_prob;

% Trajectory until t-1 or t. Preallocate for speed,
% Values are linear indices corresponding to displacement pair at each
% MC_cycle.
traj_prev = nan( (2*DY+1)* (2*DX+1), Tm*nROI );
traj_curr = traj_prev;

% Probability of the trajectories until elast completed MC_cycle. To find
% best trajectory in the end.
% P( Traj_i,j,t ) = P( del_t = i,j | Traj_i',j',t-1) * P(Traj_i',t',t-1) 
% log traj_prob( i,j,t) = log P (i,j | i',j', Traj_i', j', t-1) + log P(i',j',t-1)
traj_prob = zeros( (2*DY+1)* (2*DX+1),1);
traj_prob_new = traj_prob;

% Linear indexing of all possible displacement pairs
del_index = reshape(1:(2*DX+1)*(2*DY+1), [(2*DY+1), (2*DX+1)]);

%------------------------------------------------------------------------%
%               MOTION CORRECTION: Each Viterbi Probability
%------------------------------------------------------------------------%
tic

% Reference history for each patch upto refscale scan cycles (i.e. refscale
% * nROI) MC cycles before - average and use as new mean.

RefHist = zeros(size(Fy,2),size(Fx,2),(2*DX+1)*(2*DY+1), nROI);
tempref = zeros(size(Fy,2), size(Fx,2), (2*DX+1)*(2*DY+1));

for MC_cycle = 1:Tm*nROI
%     if mod(MC_cycle,1000)==1

%     end
   t = ceil( MC_cycle/nROI);
%    if t>refscale
%            MC_cycle
%    end
   plane = mod(MC_cycle-1, nROI) + 1;

   if t<=refscale   
    ref = R(:,:,plane);
   end
   % Fo each possible displacement( dy(t), dx(t) ):
%    tic
   for dx_t = 1:2*DX+1
   for dy_t = 1:2*DY+1
        maxlogP = -Inf;
        del_tminus1 = nan;
        
        % Reference is fixed, so one value of P1 per current displacement -
        % independent of previous displacements 
        if t<=refscale
            logP1 = sum(sum( probIgivenR( I(Fy+ del_y(dy_t), Fx+ del_x(dx_t), t, plane), ref, gamma) ) );
        end
        
        for dx_tminus1 = 1:2*DX+1
        for dy_tminus1 = 1:2*DY+1
            
            
       % Compute (log) Probability for current displacement = (i, j) for
       % each displacement pair (i', j') in previous step. Choose
       % displacement at t-1 [del(t-1)] for del(t) = (i,j) corresponding to
       % the evaluated current state that maximizes probability of
       % trajectory passing through (i,j) at time t. This is the 
       % log viterbi_traj_prob( i,j,t) 	i.e. v(i,j,t)
       
       % v(i,j,t) = log P(i, j, Traj_i,j,t)
       % v(i,j,t) = max_{i',j'}     < log P (i,j, I | i',j', Traj_i',j',t-1, R ) + log P( i', j', Traj_i',j',t-1) >
       %          = max_{i',j'}     < log P(i,j, I | i',j', R)   [MARKOV]        + v(i',j',t-1) >
       %          = max_{i',j'}     < log P( I(t)| R, i,j ) + log T(i',j'-->i,j) + v(i',j',t-1)  >
       
       % Term 1:    log P(I(t) | R, i, j) = sum_{x,y} [ log P ( ObsI(x,y,plane) | RefR(x+i, y+j, plane) ) ]
       % Term 2:    log T(i', j' --> i,j)
       % Term 3:    v(i',j',t-1)
%        

           %-----------------------------------------------------------
           % Reference is the average of last RefScale MC-corrected (for
           % current trajectory) subframes of this patch
           %-----------------------------------------------------------
    
           %%%%% Method 1: Takes too long!
%            if t==1
%                ref = R(:,:,plane);
%            else
%                ref=[];
%                tmin = max(1, t-refscale);       %Refscale previously or 1 (until enough timesteps)
%                tmax = t-1;                      % Till previous scan cycle
%                for jj=tmin:tmax
%                    dx = ceil( traj_curr( del_index(dy_tminus1, dx_tminus1), plane + (jj-1)*nROI ) / (2*DY+1) ) - DX-1 ;
%                    dy = mod( traj_curr( del_index(dy_tminus1, dx_tminus1), plane + (jj-1)*nROI ) -1,  (2*DY+1) ) - DY ;
%                    ref = cat(3, ref, squeeze(I(Fy+ dy, Fx+ dx, jj, plane)) );
%                end
%                ref = mean(ref, 3);
%            end

           %%%%% Method 2: Linear in T, in no. of displacement pair - and
           %%%%% constant for all refcale - history dependence!! increases
           %%%%% storage.
           % Keep a running sum of last k = 'rehfscale' MC-corrected intensities
%            tic
           if t>refscale
              ref = ceil(squeeze( RefHist( :,:, del_index(dy_tminus1, dx_tminus1), plane)/refscale )); 
%             logP1 = sum(sum( probIgivenR( I(Fy+ del_y(dy_t), Fx+ del_x(dx_t), t, plane), ref, gamma) ) );
              logP1 = gamma * sum(sum( squeeze(I(Fy+ del_y(dy_t), Fx+ del_x(dx_t), t, plane)).* log( ref) - ref ) );
           
           end
%            toc
           %Format of transition matrix is currently dy', dy, dx', dx
           logP2 = -r( dy_tminus1, dy_t, dx_tminus1, dx_t )/lambda; 
           
%            toc
           logP = logP1 + logP2 + traj_prob( del_index(dy_tminus1, dx_tminus1) );
           
           if logP > maxlogP
               maxlogP = logP;
               del_tminus1 = del_index(dy_tminus1, dx_tminus1);
           elseif logP == maxlogP
               del_tminus1 = [del_tminus1, del_index(dy_tminus1, dx_tminus1)];
           end
           
        end
        end
        
        %%% If multiple del(t-1) give same logP
        if size(del_tminus1,2) > 1
           % choose which one maximizes 2D-corr      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           del_tminus1 = del_tminus1(1);
        end
        
        traj_prob_new( del_index(dy_t, dx_t) ) = maxlogP;
%         try
        traj_curr( del_index(dy_t, dx_t), 1:MC_cycle ) = cat(2, traj_prev(del_tminus1, 1:MC_cycle-1), del_index(dy_t, dx_t) ); 
%         catch 
%             del_index(dy_t,dx_t)
%         end
        % Update running sum by adding best MCotion-corrected frame for
        % each trajectory
        tempref(:,:,del_index(dy_t, dx_t) ) = RefHist(:,:,del_tminus1, plane ) + ...
                                                I(Fy+dy_t-DY-1, Fx+dx_t-DX-1,t,plane);
%      toc                                     
   end
   end

   traj_prob = traj_prob_new;           % Update and set to v(i,j,t)
   traj_prev(:,1:MC_cycle) = traj_curr(:,1:MC_cycle); % Update to trajectory upto t-1

   % After computing all new trajectories, update history by removing first
   % frame and adding current(new) frame for next round of MC

   for lastdisp = 1:((2*DX+1)*(2*DY+1))
       if t>refscale
       disp1 = traj_curr( lastdisp, MC_cycle-refscale*nROI );
       dx1 = ceil(disp1/(2*DY+1) ) - DX-1;
       dy1 = mod(disp1-1, 2*DY+1) - DY;
       
           RefHist(:,:, lastdisp, plane) = tempref(:,:,lastdisp) - I(Fy+dy1, Fx+dx1,t-refscale,plane);
       else
           RefHist(:,:, lastdisp, plane) = tempref(:,:,lastdisp) ; %- I(Fy+dy1, Fx+dx1,t-refscale,plane);
       end
   end

end

%------------------------------------------------------------------------%
%           Find Best Trajectory
%------------------------------------------------------------------------%
maxlogP = max(traj_prob);
best_traj = traj_curr( (traj_prob == maxlogP), 1:Tm*nROI );

maxlogP = maxlogP - Tm*nROI*log(2*pi*lambda*lambda);  %Add the lambda-dependent contribution

toc

end