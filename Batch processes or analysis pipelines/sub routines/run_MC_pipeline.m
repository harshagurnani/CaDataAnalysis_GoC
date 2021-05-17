%% ---------------------------- SMALL PATCHES ----------------------------
if strcmpi( Basic.ROI_type, 'patch')

         r = reshape(r, [nL, nP, nTm*nTr, nROI]);
         t = reshape( t, [nL, nP, nTm*nTr, nROI] );

        if do_ref
        % No reference exists - create one 
        % -----------------------------------------------------------

                % 1. choose ROI
                % -----------------
                if Basic.MC_select_ROI
                    % Ask user to select ROI# as ref for motion correction
                    figure; nr = ceil(sqrt(nROI));
                    for roi=1:nROI,     subplot(nr,nr,roi);     imagesc( nanmean(r(:,:,:,roi), 3 ), [50 500]);     title(num2str(roi));     end
                    MC.ref_ROI  = input('Enter ROI numbers (as vector) to use for motion correction. Return empty if no need to correct:   ');
% %                     MC.ref_ROI = [1 2 3 5 6 10 28 40 11 13 15 16];
                elseif ~isfield( MC, 'ref_ROI' )
                    MC.ref_ROI = 1:nROI;
                else
                    disp(['Using ROI ' , num2str(MC.ref_ROI'), ' for MC' ] )
                end
                if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0, MC.ref_ROI = 1:nROI;  end
                if isempty(MC.ref_ROI), Basic.MC = false; end
                % -----------------------------------------------------------


                % 2. fix sess one and create ref
                % ------------------------------
                if ~isempty(MC.ref_ROI)
                    size_needed = length(MC.ref_ROI)*nL*nP*nTm*nTr*16 / (1024*1024*1024);  % gigabytes
                    if size_needed < 0.2
                        %create ref and get offsets in one shot
                        [MC.ref, MC.(sess).offsets   ] = do_mc_fastcorr( r(:,:,:,MC.ref_ROI) ); 
                        %apply offsets to all patches
                        for roi=1:nROI
                            [ ~,~, r(:,:,:,roi)]               = do_mc_fastcorr( r(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp );
                            [ ~,~, t(:,:,:,roi)]               = do_mc_fastcorr( t(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp );
                        end
                    else        % do per trial - too much data to load simultaneously
                        r = reshape(r, [nL, nP, nTm, nTr, nROI]);
                        t = reshape(t, [nL, nP, nTm, nTr, nROI]);
                        MC.(sess).offsets = cell(nTr,1);
                        % -------- create ref using trial 1
                        trial = 1;
                        [MC.ref, MC.(sess).offsets{1}   ] = do_mc_fastcorr( r(:,:,:,1,MC.ref_ROI) );
                        % ------- apply offsets on trial 1
                        for roi=1:nROI
                            [ ~,~, r(:,:,:,1,roi)]               = do_mc_fastcorr( r(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
                            [ ~,~, t(:,:,:,1,roi)]               = do_mc_fastcorr( t(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
                        end
                        % ------- get and apply offsets on other trials
                        for trial = 2:nTr
                            % get offsets
                            [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                            % apply offsets
                            for roi=1:nROI
                                [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                                [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                            end

                        end
                        r = reshape(r, [nL, nP, nTm*nTr, nROI]);
                        t=reshape( t, [nL, nP, nTm*nTr, nROI] );
                        MC.(sess).offsets = cell2mat(MC.(sess).offsets);
                    end
                end
                do_ref = false;
                % -----------------------------------------------------------


        else
        % Align to same reference
        % -----------------------------------------------------------

                if length(MC.ref_ROI)==1 && MC.ref_ROI(1)==0, MC.ref_ROI = 1:nROI;  end
                size_needed = length(MC.ref_ROI)*nL*nP*nTm*nTr*16 / (1024*1024*1024);  % gigabytes
                if ~isempty(MC.ref_ROI)
                        if size_needed < 0.2
                            % get offsets
                            [~, MC.(sess).offsets   ] = do_mc_fastcorr( r(:,:,:,MC.ref_ROI), 'ref', MC.ref );
                            % apply offsets to all patches
                            for roi=1:nROI
                                [ ~,~, r(:,:,:,roi)]                 = do_mc_fastcorr( r(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp );
                                [ ~,~, t(:,:,:,roi)]                 = do_mc_fastcorr( t(:,:,:,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets, 'max_disp', MC.max_disp );
                            end
                        else
                            r = reshape(r, [nL, nP, nTm, nTr, nROI]);
                            t = reshape(t, [nL, nP, nTm, nTr, nROI]);
                            MC.(sess).offsets = cell(nTr,1);
                            % get and apply offsets per trial
                            for trial = 1:nTr
                                %get offsets
                                [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( r(:,:,:,trial,MC.ref_ROI), 'ref', MC.ref );
                                % apply offsets
                                for roi=1:nROI
                                    [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                                    [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                                end

                            end
                            r = reshape(r, [nL, nP, nTm*nTr, nROI]);
                            t = reshape(t, [nL, nP, nTm*nTr, nROI] );
                            MC.(sess).offsets = cell2mat(MC.(sess).offsets);
                        end   
                end
        end

        % -----------------------------------------------------------------------------------------------------------------------------------------------

elseif strcmpi( Basic.ROI_type, 'plane' )
%% ----------------------------- LARGE PLANE (PATCHES) -------------------  

        if do_ref
        % No reference exists - create one 
        % -----------------------------------------------------------
            
                % draw ROI on planes for ref
                f1=figure; fixed_roi = false; this_roi = 'Yes'; ctr = 1; roi_for_ref = [];
                for roi=1:nROI
                     imagesc( nanmean(r(:,:,:,roi), 3 ));  title(num2str(roi)); hold on;
                     % get small rect patches of interest
                     while strcmpi( this_roi, 'Yes' )
                         this_roi = questdlg( 'MC Reference selection','Choose (more) ROI on current plane to use for MC ? \n (Try and select bright, prominent soma/beads)');
                         if ~strcmpi( this_roi, 'Yes' ), continue; end
                         if fixed_ROI
                             % pts as midpts for rect of fixed ROI size
                             [rx,ry] = getpts;
                             [ rect, plotrect] = make_rect( rx, ry, w/2, h/2, [nL, nP] );
                             plot( plotrect(:,1), plotrect(:,2), 'k', 'linewidth', 1 );
                             roi_for_ref = [roi_for_ref; roi, rect ];
                             ctr = ctr+1;
                         else
                             % draw rectangle
                             [rx, ry, w, h ] = getrect;
                             [ rect, plotrect] = make_rect( rx+w/2, ry+w/2, w/2, h/2, [nL, nP] );
                             plot( plotrect(:,1), plotrect(:,2), 'k', 'linewidth', 1 );
                             fixed_ROI = true;
                             roi_for_ref = [roi_for_ref; roi, rect ];
                             ctr = ctr+1;
                         end   
                     end
                     hold off;

                end
                close(f1); MC.ref_ROI = roi_for_ref;

                if ctr==1, Basic.MC = false; end
                w = MC.ref_ROI(3)-MC.ref_ROI(2)+1; h = MC.ref_ROI(5)-MC.ref_ROI(4)+1; ctr = size(MC.ref_ROI, 1);
                
                if Basic.MC

                % Crop patches
                temp_r = nan( w, h, nTm, nTr, ctr);
                for roi=1:ctr-1
                   temp_r(:,:,:,roi) = r( MC.ref_ROI(roi,2):MC.ref_ROI(roi,3), MC.ref_ROI(roi,4):MC.ref_ROI(roi,5), :, 1, MC.ref_ROI(roi,1) );
                end

                % Create ref and offsets from trial 1
                [MC.ref, MC.(sess).offsets{1}   ] = do_mc_fastcorr( temp_r(:,:,:,1,:) );
                % ------- apply offsets to large planes on trial 1
                for roi=1:nROI
                    [ ~,~, r(:,:,:,1,roi)]               = do_mc_fastcorr( r(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
                    [ ~,~, t(:,:,:,1,roi)]               = do_mc_fastcorr( t(:,:,:,1,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{1}, 'max_disp', MC.max_disp );
                end
                % ------- apply offsets to large planeson other trials
                for trial = 2:nTr
                    % get offsets
                    [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( temp_r(:,:,:,trial,:), 'ref', MC.ref );
                    % apply offsets
                    for roi=1:nROI
                        [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                        [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                    end

                end
                
                r = reshape(r, [nL, nP, nTm*nTr, nROI]);
                t=reshape( t, [nL, nP, nTm*nTr, nROI] );
                MC.(sess).offsets = cell2mat(MC.(sess).offsets);

                clear temp_r
                end
                do_ref = false;
        else

        % Align to same reference
        % -----------------------------------------------------------

            % Crop patches
            w = MC.ref_ROI(3)-MC.ref_ROI(2)+1; h = MC.ref_ROI(5)-MC.ref_ROI(4)+1; ctr = size(MC.ref_ROI, 1);
            temp_r = nan( w, h, nTm, nTr, ctr);
            for roi=1:ctr
               temp_r(:,:,:,roi) = r( MC.ref_ROI(roi,2):MC.ref_ROI(roi,3), MC.ref_ROI(roi,4):MC.ref_ROI(roi,5), :, 1, MC.ref_ROI(roi,1) );
            end

            % ------- apply offsets to large planes per trial
            for trial = 1:nTr
                % get offsets
                [~, MC.(sess).offsets{trial}   ] = do_mc_fastcorr( temp_r(:,:,:,trial,:), 'ref', MC.ref );
                % apply offsets
                for roi=1:nROI
                    [ ~,~, r(:,:,:,trial,roi)]               = do_mc_fastcorr( r(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                    [ ~,~, t(:,:,:,trial,roi)]               = do_mc_fastcorr( t(:,:,:,trial,roi), 'do_crop', true, 'use_offsets', MC.(sess).offsets{trial}, 'max_disp', MC.max_disp );
                end
            end
            
            r = reshape(r, [nL, nP, nTm*nTr, nROI]);
            t=reshape( t, [nL, nP, nTm*nTr, nROI] );
            MC.(sess).offsets = cell2mat(MC.(sess).offsets);

            clear temp_r

        end
        
        
end



if ~Basic.MC,  MC.(sess).offsets = zeros(nTm*nTr,2); end
% check for error frames
MC.(sess).error = ( abs( MC.(sess).offsets(:,1) ) > MC.max_disp | abs( MC.(sess).offsets(:,2) ) > MC.max_disp );



%% HELPER
%------------------------------------------------------------------------%

function [patch, rect_coord] = make_rect( midx, midy, halfw, halfh, roi_size )
% slice coordinates (patch), and rectangle coordinates in ccw closed
% polygon list for plotting (rect_coord)
    if midx-halfw<1, midx = halfw+1; end
    if midx+halfw>roi_size(1), midx = roi_size(1)-halfw; end
    if midy-halfh<1, midy = halfh+1; end
    if midy+halfh>roi_size(2), midy = roi_size(2)-halfh; end
    
    patch      = [ midx-halfw, midx+halfw, midy-halfh, midy+halfh ]; 
    
    rect_coord = [ midx-halfw, midy-halfh;...
                   midx+halfw, midy-halfh;...
                   midx+halfw, midy+halfh;...
                   midx-halfw, midy+halfh;...
                   midx-halfw, mody-halfh ];
end
