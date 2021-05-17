%% Get all correlations
% 
clear correlations pw_distance
for jj=1:nFiles
L = load( all_analysis{jj}, 'allAnalysed');
correlations{jj,1} = L.allAnalysed.correlations;
pw_distance{jj,1} = L.allAnalysed.pw_distance;
end

SummaryData= struct('correlations', correlations, 'pw_distance', pw_distance);

%%
allCorrelations.crus = [];
allCorrelations.lob45 = [];

for jj=find(isCrus)'
    ngoc = size(SummaryData(jj).pw_distance,1);
    usemat = triu( ones(ngoc), 1);
    allCorrelations.crus = [allCorrelations.crus; ...
                            SummaryData(jj).correlations.all.pw.corr( usemat(:)==1 ),  SummaryData(jj).pw_distance( usemat(:)==1 ), jj*ones(sum(usemat(:)),1)];
end

for jj=find(isLob45)'
    ngoc = size(SummaryData(jj).pw_distance,1);
    usemat = triu( ones(ngoc), 1);
    allCorrelations.lob45 = [allCorrelations.lob45; ...
                            SummaryData(jj).correlations.all.pw.corr( usemat(:)==1 ),  SummaryData(jj).pw_distance( usemat(:)==1 ), jj*ones(sum(usemat(:)),1)];
    
end

%% Plotting
clear a
bins = [-1:0.05:1];
[h1, bins]=hist(allCorrelations.crus(:,1),bins);
[h2, bins]=hist(allCorrelations.lob45(:,1),bins);

figure; hold on;
a(1)=area(bins,100*h1/length(allCorrelations.crus)); a(1).FaceAlpha = 0.3; a(1).FaceColor='b';
a(2)=area(bins,100*h2/length(allCorrelations.lob45)); a(2).FaceAlpha = 0.3; a(2).FaceColor='r';
xticks([-1:.5:1])
xlabel('Pairwise correlation')
ylabel('% cell pairs')
m1 = median(allCorrelations.crus(:,1));
m2 = median(allCorrelations.lob45(:,1));
plot( [m1,m1],[0,15],'b--')
plot( [m2,m2],[0,15],'m--')
legend( a, 'Crus', 'Lob IV/V')