nROI = size(allspikes,2);
coincident = zeros( nROI);

for roi1 = 1:nROI
    s = sum(allspikes(:,roi1));
    for jj = find(allspikes(:,roi1)==1)'
        tmp = allspikes(jj,:);
        try
            tmp = tmp + allspikes(jj-1,:);
        catch
            tmp = tmp + allspikes(jj+1,:);
        end
        coincident(roi1, : ) = coincident(roi1, : ) + (tmp>0);
        
        
    end
    coincident(roi1,:) = coincident(roi1,:)/s;
end