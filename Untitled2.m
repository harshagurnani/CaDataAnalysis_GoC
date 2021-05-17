noPC_Corr.all = [];
noPC_Corr.whisk = [];
for animal=1:nAnimals
for session = 1:nSess
    sess='combined';
    try
    allC    = SummaryData{animal,session}.noPC_Corr.pw.(sess).all.corr;
    nbins = floor(size(allC,3)/2);
    allC = allC(:,:,nbins+1);
    allCSig = SummaryData{animal,session}.noPC_Corr.pw.(sess).all.corr_sig(:,:,nbins+1);
    whiskC    = SummaryData{animal,session}.noPC_Corr.pw.(sess).whisk.corr;
    whiskCSig = SummaryData{animal,session}.noPC_Corr.pw.(sess).whisk.corr_sig;

    id = [(abs(allCSig)==1) & (abs(whiskCSig)==1)];
    noPC_Corr.all = [noPC_Corr.all; allC(id)];
    noPC_Corr.whisk =[noPC_Corr.whisk; whiskC(id)];
    
    
    % beh_corr
    ns = SummaryData{animal,session}.exp.nSess;
    stypes = SummaryData{animal,session}.exp.session_types;
    use_sess = SummaryData{animal,session}.exp.use_sessions;
    sessC = [];
    for ss = 1:ns
        sname = stypes{use_sess(ss)};
        allC_loco = SummaryData{animal,session}.noPC_Corr.beh.(sname).corr;
        
        
        
    catch
    end
    
end
end

%% beh corr

noPC_Corr.all = [];
noPC_Corr.whisk = [];
for animal=1:nAnimals
for session = 1:nSess
    sess='combined';
    try
    allC    = SummaryData{animal,session}.noPC_Corr.pw.(sess).all.corr;
    nbins = floor(size(allC,3)/2);
    allC = allC(:,:,nbins+1);
    allCSig = SummaryData{animal,session}.noPC_Corr.pw.(sess).all.corr_sig(:,:,nbins+1);
    whiskC    = SummaryData{animal,session}.noPC_Corr.pw.(sess).whisk.corr;
    whiskCSig = SummaryData{animal,session}.noPC_Corr.pw.(sess).whisk.corr_sig;

    id = [(abs(allCSig)==1) & (abs(whiskCSig)==1)];
    noPC_Corr.all = [noPC_Corr.all; allC(id)];
    noPC_Corr.whisk =[noPC_Corr.whisk; whiskC(id)];
    
    
    % beh_corr
    ns = SummaryData{animal,session}.exp.nSess;
    stypes = SummaryData{animal,session}.exp.session_types;
    use_sess = SummaryData{animal,session}.exp.use_sessions;
    sessC = [];
    for ss = 1:ns
        sname = stypes{use_sess(ss)};
        allC_loco = SummaryData{animal,session}.noPC_Corr.beh.(sname).corr;
        
        
        
    catch
    end
    
end
end