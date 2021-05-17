crus = 'D:\Work\OneDrive - University College London\pubs and work\Golgi in vivo imaging\Paper\Datasets\Crus';
lob45='D:\Work\OneDrive - University College London\pubs and work\Golgi in vivo imaging\Paper\Datasets\Lob4_5';
ff = dir([crus,'\*.mat']); nCrus = length(ff);
ff2 = dir([lob45,'\*.mat']); nLob45 = length(ff2);
nFiles = nCrus+nLob45;
all_analysis=cell(nFiles,1); [ isCrus, isLob45] = deal(false(nFiles,1)); 
for ctr = 1:nCrus, all_analysis{ctr} = [ff(ctr).folder,'\',ff(ctr).name]; isCrus(ctr)=true; end
for ctr = 1:nLob45, all_analysis{ctr+nCrus} = [ff2(ctr).folder,'\',ff2(ctr).name]; isLob45(ctr+nCrus) = true; end