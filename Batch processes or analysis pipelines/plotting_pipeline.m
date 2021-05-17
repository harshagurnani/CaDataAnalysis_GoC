%% PL1 - Plot Session Traces with Beh Epochs and Mean Images


for sn = 1:nSess
    sess            = session_types{use_sessions(sn)};
    nROI = size(Norm.(sess).n,2);

    if contains(sess, 'AP') || contains(sess,'puff') && get_puff
        add_puff = true;
    else, 
    add_puff = false; 
    end
    
%     % Traces and Patches only
%     if get_plot_simple
%         show_patch_and_trace
%         savefig([exp.analysed,'/',sess,'_norm.fig']);
%     end
    
%     With Beh Periods and MI
    if get_plot_beh
        show_patch_and_videotrace
        savefig([exp.analysed,'/',sess,'_vidMI.fig']);
    end
    
end

clear cj plot_exp 
%------------------------------------------------------------------------------------------------------------------------------------

%% Plot cells in 3D
for sn = 1
   sess = session_types{use_sessions(sn)};
   plot_patch_in_volume(p.POI_coordinates, 'z');
   savefig([exp.analysed,'/',sess,'_patch_in_vol.fig']);

end

%% Plot masked mean images
masks=Seg.masks; nROI = numel(masks);
nr = ceil(sqrt(nROI));

% figure;for jj=1:nROI, subplot(nr,nr,jj); imagesc(masks{zid(jj)}); end
% savefig([exp.analysed,'\Masks.fig'])
figure;colormap(gray);for jj=1:nROI, subplot(nr,nr,jj); imagesc(Seg.mean_image{jj}); yticks([]); xticks([]); end
savefig([exp.analysed,'/Mean image.fig'])
figure;colormap(gray);for jj=1:nROI, subplot(nr,nr,jj); imagesc(Seg.mean_image{jj}.*Seg.masks{jj}); xticks([]); yticks([]); end
savefig([exp.analysed,'/soma masks.fig'])
figure;colormap(gray);for jj=1:nROI, subplot(nr,nr,jj); imagesc(Seg.mean_image{jj}.*Seg.background_masks{jj}); xticks([]); yticks([]); end
savefig([exp.analysed,'/dark background masks.fig'])
    


%% Backgrounds

% Background masks
% figure;colormap(gray); nr=ceil(sqrt(nROI)); for jj=1:nROI, subplot(nr,nr,jj); imagesc(Seg.background_masks{jj}.*Seg.mean_image{jj});yticks([]); xticks([]);end

[~,zid] = sort( p.POI_coordinates.ROIs_Z(:,1));
for sn =1:nSess
    sess = session_types{use_sessions(sn)};
    xx = Norm.(sess).background.*repmat(Norm.(sess).background_F0./Norm.(sess).F0, p.nTimepoints*p.nTrials, 1);
    
    figure; colormap(jet);
    plot(Norm.(sess).t(:,zid)/1000, xx(:,zid)+(1:nROI)*2,           'w--','linewidth', 1.2); hold on; 
    plot(Norm.(sess).t(:,zid)/1000, Norm.(sess).n(:,zid)+(1:nROI)*2,'linewidth', 1.2); 
    xlabel('Time (s)');
    suptitle(['Session ' sess]);
    savefig([exp.analysed, '\Soma and background_', sess,'.fig'])
end

%% Correlation Plots

% Plot lagcorr of individual soma, and somatic components

% C with Whisking
sess = 'baseline'; beh = 'whisk';
do_lagcorr_plots_all;


sess = 'baseline'; beh = 'loco';
do_lagcorr_plots_all;


if exist('Pupil', 'var')
sess = 'baseline'; beh = 'pupil';
do_lagcorr_plots_all;


end

if isfield(Corr,'AP')
    sess = 'AP'; beh = 'whisk';
    do_lagcorr_plots_all;

    sess = 'AP'; beh = 'loco';
    do_lagcorr_plots_all;
end


disp('.....................Saved Plots')

%% Plot Responses to Air Puff
plot_AP_responses2
