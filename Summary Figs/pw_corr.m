%% pw corr
Corr.all = [];
Corr.whisk = [];
BehCorr.loco = [];
BehCorr.whisk = [];
BehCorr.Id = [];
BehCorr.zeroloco = [];
BehCorr.zerowhisk = [];

for animal=1:nAnimals
for session = 1:nSess
    sess='combined';
    if ~isempty(SummaryData{animal,session})
%     try
    allC    = SummaryData{animal,session}.Corr.pw.(sess).all.corr;
    nbins = floor(size(allC,3)/2);
    allC = allC(:,:,nbins+1);
    allCSig = SummaryData{animal,session}.Corr.pw.(sess).all.corr_sig(:,:,nbins+1);
    whiskC    = SummaryData{animal,session}.Corr.pw.(sess).whisk.corr;
    whiskCSig = SummaryData{animal,session}.Corr.pw.(sess).whisk.corr_sig;

    id = [(abs(allCSig)==1) & (abs(whiskCSig)==1)];
    Corr.all = [Corr.all; allC(id)];
    Corr.whisk =[Corr.whisk; whiskC(id)];
%     catch
%     end
%     
    % beh_corr
%     try    
%     ns = SummaryData{animal,session}.exp.nSess;
%     stypes = SummaryData{animal,session}.exp.session_types;
%     use_sess = SummaryData{animal,session}.exp.use_sessions;
%     sessC = [];
%     for ss = 1:ns
%         sname = stypes{use_sess(ss)};
%         allC_loco = SummaryData{animal,session}.Corr.beh.(sname).corr;
%         
%     end
    allC_loco = SummaryData{animal,session}.Corr.beh.noPC1.loco.corr;
    allCSig_loco = SummaryData{animal,session}.Corr.beh.noPC1.loco.corr_sig;
    allCloco_lags = SummaryData{animal,session}.Corr.beh.noPC1.loco.lags;
    allC_whisk = SummaryData{animal,session}.Corr.beh.noPC1.whisk.corr;
    allCSig_whisk = SummaryData{animal,session}.Corr.beh.noPC1.whisk.corr_sig;
    allCwhisk_lags = SummaryData{animal,session}.Corr.beh.noPC1.whisk.lags;
    nroi = size(allC_whisk,1);
    nbins = floor(size(allC_whisk,2)/2);
    for roi=1:nroi
        locoC = [NaN NaN];zeroloco = NaN;
         if max(abs(smooth(SummaryData{animal,session}.Speed(:,2),50)))*1.2>10 
        [pk1, id1]=findpeaks(allC_loco(roi,:),'NPeaks',1);
        [pk2, id2]=findpeaks(-allC_loco(roi,:),'NPeaks',1);
        if ~isempty(id1) && ~isempty(id2)
            if allCSig_loco(roi,id1)==1 && allCSig_loco(roi,id2)==-1 && pk1>pk2,locoC = [allC_loco(roi,id1),allCloco_lags(id1)];
            elseif allCSig_loco(roi,id1)==1 && allCSig_loco(roi,id2)==-1 && pk1<pk2,locoC = [allC_loco(roi,id2),allCloco_lags(id2)];
            elseif allCSig_loco(roi,id1)==1, locoC = [allC_loco(roi,id1),allCloco_lags(id1)];
            elseif allCSig_loco(roi,id2)==-1, locoC = [allC_loco(roi,id2),allCloco_lags(id2)];    
            end
        elseif ~isempty(id1) && allCSig_loco(roi,id1)==1, locoC = [allC_loco(roi,id1),allCloco_lags(id1)];
        elseif ~isempty(id2) && allCSig_loco(roi,id2)==-1, locoC = [allC_loco(roi,id2),allCloco_lags(id2)];
        end 
         if abs(allCSig_loco(roi,nbins+1))==1, zeroloco = allC_loco(roi,nbins+1); end
         end
        whiskC = [NaN NaN]; zerowhisk = NaN;
        [pk1, id1]=findpeaks(allC_whisk(roi,:),'NPeaks',1);
        [pk2, id2]=findpeaks(-allC_whisk(roi,:),'NPeaks',1);
        if ~isempty(id1) && ~isempty(id2)
            if allCSig_whisk(roi,id1)==1 && allCSig_whisk(roi,id2)==-1 && pk1>pk2,whiskC = [allC_whisk(roi,id1),allCwhisk_lags(id1)];
            elseif allCSig_whisk(roi,id1)==1 && allCSig_whisk(roi,id2)==-1 && pk1<pk2,whiskC = [allC_whisk(roi,id2),allCwhisk_lags(id2)];
            elseif allCSig_whisk(roi,id1)==1, whiskC = [allC_whisk(roi,id1),allCwhisk_lags(id1)];
            elseif allCSig_whisk(roi,id2)==-1, whiskC = [allC_whisk(roi,id2),allCwhisk_lags(id2)];    
            end
        elseif ~isempty(id1) && allCSig_whisk(roi,id1)==1, whiskC = [allC_whisk(roi,id1),allCwhisk_lags(id1)];
        elseif ~isempty(id2) && allCSig_whisk(roi,id2)==-1, whiskC = [allC_whisk(roi,id2),allCwhisk_lags(id2)];
        end 
       
        if abs(allCSig_whisk(roi,nbins+1))==1, zerowhisk = allC_whisk(roi,nbins+1); end
        
        BehCorr.loco = [BehCorr.loco; locoC];
        BehCorr.whisk = [BehCorr.whisk; whiskC];
        BehCorr.Id = [BehCorr.Id; animal, session, roi];
        BehCorr.zeroloco = [BehCorr.zeroloco; zeroloco];
        BehCorr.zerowhisk = [BehCorr.zerowhisk; zerowhisk];
    end
%     catch
%     end
    end
end
end

