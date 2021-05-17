allData.neurons.f = activity_all.soma;
allData.neurons.time = activity_all.t;
allData.ROI.keep_f = 1:p.nROIs;
allData.ROI.keep_spikes = gid;
allData.neurons.spikes = activity_all.spikes;
dt = 1/p.acquisition_rate;
allspikes = zeros(size(allData.neurons.spikes));
for jj=1:length(spk_manual), for kk=1:length(spk_manual(jj).spikes), tmp = ceil(spk_manual(jj).spikes(kk)/dt); 
        allspikes(tmp,jj)=allspikes(tmp, jj)+1; end
end
allData.neurons.spikes = allspikes;