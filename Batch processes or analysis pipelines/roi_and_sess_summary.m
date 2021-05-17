


for bn = 1:numel(allbeh)
% sesstype = 'baseline'; beh = 'whisk';
beh = allbeh{bn};
allROI_behcorr.(beh) = [];
allsess_corr.(beh)=[];
allROI_behcorrlag.(beh)=[];

    for animal = 1:nAnimals
        nSess = numel(SummaryData.BehCorr{animal});
        for session = 1:nSess
            checkon=true;
            if ~isempty( SummaryData.BehCorr{animal}{session} )
                
                nROI = length(StatData.BehCorr{animal}{session}.baseline.(beh).maxC);
                sesstemp = [];
                for roi=1:nROI
                    %Baseline
                    tmp = []; tmp_all =[]; 
                    if isfield(StatData.BehCorr{animal}{session},'baseline')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.baseline.(beh).maxC(roi)];
                    lagtmp=StatData.BehCorr{animal}{session}.baseline.(beh).lagC(roi);
%                     if lagtmp<-2000
%                         lagtmp=lagtmp+2000;
%                     elseif lagtmp>2000,
%                         lagtmp = lagtmp-2000;
%                     end
                    allROI_behcorrlag.(beh)=[allROI_behcorrlag.(beh); lagtmp]; 
                    end
                    if isfield(StatData.BehCorr{animal}{session},'baseline2')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.baseline2.(beh).maxC(roi)];
                    end
                    if isfield(StatData.BehCorr{animal}{session},'baseline3')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.baseline3.(beh).maxC(roi)];
                    end
                    tmp_baseline = nanmean(tmp);
                    tmp_all=[tmp_all, tmp];
                    
                    %AP
                    tmp = []; 
                    if isfield(StatData.BehCorr{animal}{session},'AP')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.AP.(beh).maxC(roi)];
                    end
                    if isfield(StatData.BehCorr{animal}{session},'AP2')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.AP2.(beh).maxC(roi)];
                    end
                    tmp_AP = nanmean(tmp);
                    tmp_all=[tmp_all, tmp];
                    
                    %puff_other
                    tmp = []; 
                    if isfield(StatData.BehCorr{animal}{session},'puff_paw')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.AP.(beh).maxC(roi)];
                    end
                    if isfield(StatData.BehCorr{animal}{session},'puff_out')
                    tmp = [tmp, StatData.BehCorr{animal}{session}.AP2.(beh).maxC(roi)];
                    end
                    tmp_op = nanmean(tmp);
                    tmp_all=[tmp_all, tmp];
                    
                    %any
                    tmp_any = nanmean(tmp_all);
                    
                    allROI_behcorr.(beh) = [allROI_behcorr.(beh); tmp_baseline, tmp_AP, tmp_op, tmp_any];
                    sesstemp = [sesstemp;  tmp_baseline, tmp_AP, tmp_op, tmp_any];
                end
                allsess_corr.(beh)   = [allsess_corr.(beh); nanmean( sesstemp)];
                

            end
        end
    end
end
