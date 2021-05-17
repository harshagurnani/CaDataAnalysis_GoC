% create summary data

nAnimals = size(all_analysis,1); nExp = size(all_analysis,2);
nSess = size(all_analysis,3);

for an=1:nAnimals
    set = 1;
for ex = 1:nExp
for sn = 1:nSess
   
if ~isempty(all_analysis{ an, ex, sn })
    fpath = all_analysis{ an, ex, sn };
    L  = load( [ fpath, '\exp_details.mat'],'exp');
    SummaryData{an,set}.exp = L.exp;
%     L  = load( [ fpath, '\correlations.mat'],'Corr');
%     SummaryData{an,set}.Corr = L.Corr;
    
    L  = load( [ fpath, '\masked_norm_test.mat'],'Norm');
    SummaryData{an,set}.Norm = L.Norm;
    L  = load( [ fpath, '\masked_norm_test.mat'],'p');
    SummaryData{an,set}.p = L.p;
%     L  = load( [ fpath, '\PCA_results.mat'],'PCA');
%     SummaryData{an,set}.PCA = L.PCA;
%     L  = load( [ fpath, '\PC_Corr.mat'],'PC_Corr');
%     SummaryData{an,set}.PC_Corr = L.PC_Corr;
    L  = load( [ fpath, '\ROI.mat'],'ROI');
    SummaryData{an,set}.ROI = L.ROI;
%     L  = load( [ fpath, '\PuffResponses.mat'],'PCA');
%     SummaryData{an,set}.PuffResponses = L.PuffResponses;
%     L  = load( [ fpath, '\conactenated.mat'],'beh_all');
%     SummaryData{an,set}.Speed = L.beh_all.loco;
%     L  = load( [ fpath, '\conactenated.mat'],'activity_all');
%     SummaryData{an,set}.norm_soma = L.activity_all.soma;
%     L  = load( [ fpath, '\SomaResponses_noPC1.mat'],'SomaResponses');
%     SummaryData{an,set}.noPC1_resp = L.SomaResponses;
%     L  = load( [ fpath, '\corr_analysis.mat'],'CorrAnalysis');
%     SummaryData{an,set}.pwCA = L.CorrAnalysis; 
    
    set = set+1;
    
    
end
    
end
end
end
nSess = size(SummaryData,2);