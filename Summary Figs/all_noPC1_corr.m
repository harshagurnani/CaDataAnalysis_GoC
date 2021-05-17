% get correlation histogram

nAnimals = size(all_analysis,1);
nExp = size(all_analysis,2);
nSess = size(all_analysis,3);


noPC1_corr = [];


qual_zlimit = [-3 4];

for animal = 1:nAnimals
for session=1:nSess
    if ~isempty(  SummaryData{animal, session} )
       Corr = SummaryData{animal, session}.Corr ;
       Norm = SummaryData{animal, session}.norm_soma;
       
       sname = 'noPC1';
       temp = Corr.pw.(sname).all.corr;
       lags = floor(size(temp,3)/2); nROI = size(temp,1);
       temp = temp(:,:,lags+1);     
       sig = Corr.pw.(sname).all.corr_sig(:,:,lags+1);
       noresp = ( ( max(zscore(Norm))<qual_zlimit(2)) & (min(zscore(Norm))> qual_zlimit(1)) );
       %no diag, no double counting
       ind = triu( ones(nROI)) - eye(nROI);
       %remove non-resp
       ind(noresp,:) = 0; ind(:, noresp) = 0;
       temp = temp(ind==1);
       sig=sig(ind==1);
       %save corr and sig
       
       noPC1_corr = [ noPC1_corr; temp sig ];
               
           
       end
    end

end
