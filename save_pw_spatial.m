for an=1:2
allPW{an}.D.pw     =[];
allPW{an}.D.pwP    =[];
allPW{an}.D.z     = [];
allPW{an}.D.x     =[];
allPW{an}.D.y     =[];
allPW{an}.peakCC  =[];
allPW{an}.peakW   =[];
allPW{an}.peakSig = [];

allPW{an}.peakCC_n1  =[];
allPW{an}.peakW_n1   =[];
allPW{an}.peakSig_n1 = [];

allPW{an}.set     = [];


for set = 1:size(SummaryData,2)
    if ~isempty( SummaryData{an,set})
    
        % distances
        pwCA = SummaryData{an,set}.pwCA;
        allPW{an}.D.x = [ allPW{an}.D.x; pwCA.px2um*pwCA.allD.rAP];
        allPW{an}.D.y = [ allPW{an}.D.y; pwCA.px2um*pwCA.allD.rML];
        
        % peak cc
        allPW{an}.peakCC = [ allPW{an}.peakCC; pwCA.peakCC];
        allPW{an}.peakW = [ allPW{an}.peakW; pwCA.peakW];
        allPW{an}.peakSig = [ allPW{an}.peakSig; pwCA.peakSig];
        
        allPW{an}.peakCC_n1 = [ allPW{an}.peakCC_n1; pwCA.peakCC_noPC1];
        allPW{an}.peakW_n1 = [ allPW{an}.peakW_n1; pwCA.peakW_noPC1];
        allPW{an}.peakSig_n1 = [ allPW{an}.peakSig_n1; pwCA.peakSig_noPC1];
        
        allPW{an}.set    = [ allPW{an}.set; set*ones(size(pwCA.peakCC))];
        
        dz = [];
        nroi = size(pwCA.Locations.ROIs_X,1);
        for jj=1:nroi, for kk=jj+1:nroi
            dz = [dz; abs(pwCA.Locations.ROIs_Z(jj,1)-pwCA.Locations.ROIs_Z(kk,1))];
        end
        end
        allPW{an}.D.z = [allPW{an}.D.z; dz ];
    end
    
end
    
     allPW{an}.D.pwP = [allPW{an}.D.pwP; sqrt( allPW{an}.D.x.^2+ allPW{an}.D.y.^2 )];
     allPW{an}.D.pw = [allPW{an}.D.pw; sqrt( allPW{an}.D.x.^2+ allPW{an}.D.y.^2 + allPW{an}.D.z.^2 )];
end
