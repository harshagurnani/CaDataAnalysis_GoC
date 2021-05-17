%load raw
allsess = dir([exppath,'\raw_*.mat']);
for jj=1:length(allsess),
    load( [exppath, '\', allsess(jj).name]);
end