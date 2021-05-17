% given summary data


nAnimals = 2;

allsess = { 'baseline', 'AP', 'AP2', 'baseline2', 'baseline3', 'puff_out', 'puff_paw'};
allbeh  = { 'whisk', 'loco', 'pupil'};



 %% Mean behavioural correlations
for sn =1:numel(allsess)
sesstype=allsess{sn};
for bn = 1:numel(allbeh)
% sesstype = 'baseline'; beh = 'whisk';
beh = allbeh{bn};


   
    MeanCorr.(sesstype).(beh) = [];

    for animal = 1:nAnimals
        nSess = numel(SummaryData.BehCorr{animal});

        for session = 1:nSess
            
            if ~isempty( SummaryData.BehCorr{animal}{session} )
                if isfield(SummaryData.BehCorr{animal}{session}, sesstype)
                    allcorrs = SummaryData.BehCorr{animal}{session}.(sesstype);

                    % whisking
                    nROI = size(allcorrs.whisk.corr,1);
                    StatData.BehCorr{animal}{session}.(sesstype).(beh).maxC = nan(nROI,1);
                    StatData.BehCorr{animal}{session}.(sesstype).(beh).lagC = nan(nROI,1);
                    if isfield(allcorrs, beh)
                        for roi=1:nROI
                            validc  = allcorrs.(beh).corr(roi, allcorrs.(beh).corr_sig(roi,:)~=0);
                            validlag = allcorrs.(beh).lags(allcorrs.(beh).corr_sig(roi,:)~=0);

                            if ~isempty(validc)
                                [peakposc, indpos] = findpeaks( validc, validlag, 'SortStr','descend' );
                                [peaknegc, indneg] = findpeaks( -validc,validlag, 'SortStr','descend' );
                                if ~isempty(peakposc),   peakposc = peakposc(1); indpos = indpos(1); else, peakposc = 0; indpos = NaN; end
                                if ~isempty(peaknegc),   peaknegc = peaknegc(1); indneg = indneg(1); else, peaknegc = 0; indneg = NaN; end

                                if abs(peakposc)>=abs(peaknegc)
                                    peakC = peakposc; lagC = indpos;
                                else
                                    peakC = peaknegc; lagC = indneg;
                                end
                                StatData.BehCorr{animal}{session}.(sesstype).(beh).maxC(roi) = peakC;
                                StatData.BehCorr{animal}{session}.(sesstype).(beh).lagC(roi) = lagC;
                            end
                        end

                        StatData.BehCorr{animal}{session}.(sesstype).(beh).meanC = nanmean(StatData.BehCorr{animal}{session}.(sesstype).(beh).maxC);
                        MeanCorr.(sesstype).(beh) = [ MeanCorr.(sesstype).(beh); nanmean(StatData.BehCorr{animal}{session}.(sesstype).(beh).maxC)];
                    end
                end
            end
        end
    end

end
end

%% Mean PC-behavioural correlations
for sn =1:numel(allsess)
sesstype=allsess{sn};
for bn = 1:numel(allbeh)
% sesstype = 'baseline'; beh = 'whisk';
beh = allbeh{bn};


    
    TopPCCorr.(sesstype).(beh) = [];

    for animal = 1:nAnimals
        nSess = numel(SummaryData.BehCorr_PC{animal});

        for session = 1:nSess
            StatData.BehCorr_PC{animal}{session} = [];
            if ~isempty( SummaryData.BehCorr_PC{animal}{session} )
                if isfield(SummaryData.BehCorr_PC{animal}{session}, sesstype)
                    allcorrs = SummaryData.BehCorr_PC{animal}{session}.(sesstype);

                    % whisking
                    nPC = size(allcorrs.whisk.corr,1);
                    StatData.BehCorr_PC{animal}{session}.(sesstype).(beh).maxC = nan(nPC,1);
                    StatData.BehCorr_PC{animal}{session}.(sesstype).(beh).lagC = nan(nPC,1);
                    if isfield(allcorrs, beh)
                        for roi=1:nPC
                            validc  = allcorrs.(beh).corr(roi, allcorrs.(beh).corr_sig(roi,:)~=0);
                            validlag = allcorrs.(beh).lags(allcorrs.(beh).corr_sig(roi,:)~=0);

                            if length(validc)>2
                                [peakposc, indpos] = findpeaks( validc, validlag, 'SortStr','descend' );
                                [peaknegc, indneg] = findpeaks( -validc,validlag, 'SortStr','descend' );
                                if ~isempty(peakposc),   peakposc = peakposc(1); indpos = indpos(1); else, peakposc = 0; indpos = NaN; end
                                if ~isempty(peaknegc),   peaknegc = peaknegc(1); indneg = indneg(1); else, peaknegc = 0; indneg = NaN; end

                                if abs(peakposc)>=abs(peaknegc)
                                    peakC = peakposc; lagC = indpos;
                                else
                                    peakC = peaknegc; lagC = indneg;
                                end
                                StatData.BehCorr_PC{animal}{session}.(sesstype).(beh).maxC(roi) = peakC;
                                StatData.BehCorr_PC{animal}{session}.(sesstype).(beh).lagC(roi) = lagC;
                            end
                        end

                        ntop = min(5, nPC);
                        toppccorr = nan(1,5);
                        toppccorr(1:ntop) = reshape(StatData.BehCorr_PC{animal}{session}.(sesstype).(beh).maxC(1:ntop),[1,ntop]);
                        TopPCCorr.(sesstype).(beh) = [ TopPCCorr.(sesstype).(beh); toppccorr ];
                    else
                        TopPCCorr.(sesstype).(beh) = [ TopPCCorr.(sesstype).(beh); nan(1,5) ];
                    end
                end
            end
        end
    end

end
end
