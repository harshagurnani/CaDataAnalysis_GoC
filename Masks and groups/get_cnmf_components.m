function [compMatrix, binary_masks,ROI, params] = get_cnmf_components( raw, K, show_ini, select_comp, nCycles, varargin )
%% Get spatial components from patch data
% Using CalmAn Toolbox : https://github.com/flatironinstitute/CaImAn-MATLAB
% Adapted from demo_script.m in CalmAn Toolbox
%
% Inputs:
%   raw         4-D array   Raw imaging data ([x y t ROI])
%       If ROI = 1, data from single patch, else different patches
%   concat      Boolean     Concatenate patches for cnmf? (default=false)
%       If data is from small independent patches, you can concatenate by
%       adding few extra padding in between.
%   varargin    Options for CNMF
%       Can be structure created with CNMFSetParams.m, or name-value pairs
%
% Outputs:
%   masks       Cell array  Weighted masks for each pixel
%   ROI         Array       Patch(ROI) number corresponding to each mask
%
%
% Modified by: Harsha Gurnani. Dec 2018

% Create Params
[d1, d2, nT, nROI] = size( raw );
if isempty(K), K = 10; end
if isempty(show_ini), show_ini = false; end
if isempty(select_comp), select_comp = false; end
if isempty(nCycles), nCycles = 2; end
if ~isempty(varargin)
    if ~isstruct(varargin{1})
        params = CNMFSetParams( 'd1', d1, 'd2', d2, varargin{:} );
    else
        params = varargin{1}; params.d1 = d1; params.d2 = d2; params.d = d1*d2; 
    end
else
    params = make_default_params(d1,d2);
end

p = params.p;
tau = params.gSig;
d = d1*d2;


binary_masks=[];
for roi=1:nROI

    fprintf('Processing ROI %d ...\n', roi)
    Y = squeeze(raw(:,:,:,roi));
    
%      if ~isa(Y,'single');    Y = single(Y);  end                            % convert to single
    [P,Y] = preprocess_data(Y,p);                                           % Data pre-processing: Noise power spectrum, time constants and saturation
    [Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,params,P);    % initialize
   
    Cn =  correlation_image(Y);
    if show_ini %show initial guess
    figini=figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title(sprintf('Center of components in Patch %d found from initialization algorithm',roi));
    drawnow;disp('Close figure to continue...');
    waitfor(figini);
    end

    %% First CNMF
    
    %Update Spatial-1
    Yr = reshape(Y,d,nT);
    [A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,params);
    %Update Temporal-1
    P.p = 0;    % set AR temporarily to zero for speed
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,params);
    
    %% Subsequent cycles of CNMF
    nIter = 1;
    while nIter < nCycles
    
        %% Refine components

        % Correlate components
        rval_space = classify_comp_corr(Y,A,C,b,f,params);
        ind_corr = rval_space > params.space_thresh;           % components that pass the correlation test

        % Trace 'probability' given noise distribution 
    %     fitness = compute_event_exceptionality(C+YrA,params.N_samples_exc,params.robust_std);
    %     ind_exc = (fitness < params.min_fitness);   %the lesser, the better. spiking cells will have lower value of fitness

    %     keep = ind_corr & ind_exc;
        keep = ind_corr;

        %Merge components
        A_keep = A(:,keep); C_keep = C(keep,:);
        [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A_keep,b,C_keep,f,P,S,params);

        %% CNMF
        if ~isempty(A_keep)
        Pm.p = p;    % restore AR value
        %Update Spatial-2
        [A,b,C] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,params);
        %Update temporal-2
        [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,C,f,Pm,params);
        end
        nIter = nIter+1;
    end
    
    if ~isempty(A_keep)
    [A_or,C_or,S_or,P_or] = order_ROIs(A,C,S,P); % order components
    K_m = size(C_or,1); %No. of components
    else
        K_m = 0;
    end
    
    
    %% Plotting
    if K_m
    if select_comp
        figcomp{roi} = figure;
        [Coor{roi},~] = plot_contours(A_or,Cn,params,1); % contour plot of spatial footprints
        title(sprintf('Patch %d',roi));
        figcompsep{roi} = figure;
        patches=extract_patch(A_or,[d1,d2],[30,30]);
        nr = ceil(sqrt(K_m));
        for jj=1:K_m
           subplot( nr,nr,jj); imagesc(patches(:,:,1,jj));xticks([]);yticks([]);title(num2str(jj));
        end
        suptitle(sprintf('Components from Patch %d',roi));
        keep_comp = input( 'Enter component no.s to keep:');
        if keep_comp == 0, keep_comp=1:K_m; end
        close all
    else
        keep_comp=1:K_m;
    end
    nKeep = length(keep_comp);
    [pixID, CompID, Wts] = find(A_or);
    delrows = ~ismember(CompID, keep_comp);
    pixID( delrows ) = []; CompID( delrows ) = []; Wts( delrows ) = [];
    tmpcomp = CompID;
    for jj=1:nKeep, tmpcomp(CompID==keep_comp(jj)) = jj; end 
    CompID = tmpcomp;
    tmp_masks = cell(nKeep,1);
    
    for jj=1:nKeep
    mid = (CompID==jj); nm = sum(mid);
    tmp_masks{jj} = sparse(pixID( mid ), ones( nm,1 ), Wts(mid), d, 1 );
    tmp_masks{jj} = full( tmp_masks{jj});
    tmp_masks{jj} = reshape( tmp_masks{jj}, d1, d2 );
    end
    tmp_ROI=ones(nKeep,1)*roi;   
    compMatrix{roi} = sparse( pixID, CompID, Wts, d, nKeep);
    if roi==1
       binary_masks = tmp_masks;
       ROI = tmp_ROI;
    else
        binary_masks = [binary_masks; tmp_masks];
        ROI = [ROI; tmp_ROI];
    end
    end
end





function pm = make_default_params(d1,d2)
    K = 30;                                           % number of components to be found
    tau = 8;                                          % std of gaussian kernel (half size of neuron) 
    p = 2;                                            % order of AR dynamics 
    merge_threshold = 0.85;                           % merging threshold
    nBG_comp = 1;                                     % number of background components 
    min_SNR = 3;                                      % minimum SNR threshold
    scorr_thr = 0.65;                                  % space correlation threshold
    ellipse_size = [3 10];
    
    pm = CNMFSetParms(...   
        'd1',d1,'d2',d2,...                         % dimensionality of the FOV
        'p',p,...                                      
        'gSig',tau,...                              
        'merge_thr',merge_threshold,...                         
        'nb',nBG_comp,...                                     
        'min_SNR',min_SNR,...
        'space_thresh',scorr_thr,...
        'cnn_thr',0.2,...                            % threshold for CNN classifier    
        'min_size', ellipse_size(1),...
        'max_size', ellipse_size(2) ...
        );
