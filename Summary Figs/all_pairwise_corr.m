% get correlation histogram

nAnimals = size(all_analysis,1);
nExp = size(all_analysis,2);
nSess = size(all_analysis,3);

all_baseline_corr = [];
all_combined_corr = [];
all_AP_corr = [];
all_noPC1_corr = [];

max_baseline_corr = [];
max_combined_corr = [];
max_AP_corr = [];


qual_zlimit = [-3 4];

for animal = 1:nAnimals
for session=1:nSess
    if ~isempty(  SummaryData{animal, session} )
       Corr = SummaryData{animal, session}.Corr ;
       Norm = SummaryData{animal, session}.Norm;
       Norm = SummaryData{animal, session}.Norm;
       p = SummaryData{animal,session}.p;
       ns = fieldnames( Corr.pw);
       for stype = 1:length(ns)
           sname = ns{stype};
           if contains( sname, 'baseline')
               temp = Corr.pw.(sname).all.corr;
               [maxcorr, mincorr] = deal(nan(p.nROIs));
               lags = floor(size(temp,3)/2); nROI = size(temp,1);
               noresp = ( ( max(zscore(Norm.baseline.n))<qual_zlimit(2)) & (min(zscore(Norm.baseline.n))> qual_zlimit(1)) );
               ind = triu( ones(nROI)) - eye(nROI);
               ind(noresp,:) = 0; ind(:, noresp) = 0;
               maxcorr = max(temp.*(Corr.pw.(sname).all.corr_sig==1),3);
               mincorr = min(temp.*(Corr.pw.(sname).all.corr_sig==-1),3);
               maxcorr(abs(mincorr)> maxcorr)= mincorr(abs(mincorr)> maxcorr);
               temp = temp(:,:,lags+1);              
               maxcorr = maxcorr(ind == 1 );
               ind(Corr.pw.(sname).all.corr_sig(:,:,lags+1)==0)=0;
               temp = temp(ind == 1);
               all_baseline_corr = [ all_baseline_corr; temp ];
               max_baseline_corr = [ max_baseline_corr; temp ];

           elseif contains( sname, 'AP')
               
               temp = Corr.pw.(sname).all.corr;
               [maxcorr, mincorr] = deal(nan(p.nROIs));
               lags = floor(size(temp,3)/2); nROI = size(temp,1);
               noresp = ( ( max(zscore(Norm.baseline.n))<qual_zlimit(2)) & (min(zscore(Norm.baseline.n))> qual_zlimit(1)) );
               ind = triu( ones(nROI)) - eye(nROI);
               ind(noresp,:) = 0; ind(:, noresp) = 0;
               maxcorr = max(temp.*(Corr.pw.(sname).all.corr_sig==1),3);
               mincorr = min(temp.*(Corr.pw.(sname).all.corr_sig==-1),3);
               maxcorr(abs(mincorr)> maxcorr)= mincorr(abs(mincorr)> maxcorr);
               temp = temp(:,:,lags+1);              
               maxcorr = maxcorr(ind == 1 );
               ind(Corr.pw.(sname).all.corr_sig(:,:,lags+1)==0)=0;
               temp = temp(ind == 1);
               all_AP_corr = [ all_AP_corr; temp ];
               max_AP_corr = [ all_AP_corr; temp ];
               
           elseif contains( sname, 'combined')
               
               temp = Corr.pw.(sname).all.corr;
               [maxcorr, mincorr] = deal(nan(p.nROIs));
               lags = floor(size(temp,3)/2); nROI = size(temp,1);
               noresp = ( ( max(zscore(Norm.baseline.n))<qual_zlimit(2)) & (min(zscore(Norm.baseline.n))> qual_zlimit(1)) );
               ind = triu( ones(nROI)) - eye(nROI);
               ind(noresp,:) = 0; ind(:, noresp) = 0;
               maxcorr = max(temp.*(Corr.pw.(sname).all.corr_sig==1),3);
               mincorr = min(temp.*(Corr.pw.(sname).all.corr_sig==-1),3);
               maxcorr(abs(mincorr)> maxcorr)= mincorr(abs(mincorr)> maxcorr);
               temp = temp(:,:,lags+1);              
               maxcorr = maxcorr(ind == 1 );
%                ind(Corr.pw.(sname).all.corr_sig(:,:,lags+1)==0)=0;
               temp = temp(ind == 1);
               all_combined_corr = [ all_combined_corr; temp ];
               max_combined_corr = [ max_combined_corr; temp ];
               
           elseif contains( sname, 'noPC1')
               
               temp = Corr.pw.(sname).all.corr;
               lags = floor(size(temp,3)/2); nROI = size(temp,1);
               temp = temp(:,:,lags+1);     
               
               noresp = ( ( max(zscore(Norm.baseline.n))<qual_zlimit(2)) & (min(zscore(Norm.baseline.n))> qual_zlimit(1)) );
               %no diag, no double counting
               ind = triu( ones(nROI)) - eye(nROI);
               %remove non-resp
               ind(noresp,:) = 0; ind(:, noresp) = 0;
%                %remove non-dig
%                ind(Corr.pw.(sname).all.corr_sig(:,:,lags+1)==0)=0;
               temp = temp(ind == 1);
               all_noPC1_corr = [ all_noPC1_corr; temp ];
               
           end
       end
    end

end
end