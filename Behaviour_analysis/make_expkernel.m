%% make kernel

tau1=40;
t0 = [-100 0 100];
tau2 = [40 100 150];
% x = -200:20:6000;
dt = 1000/300;
x = -30*dt:dt:1000;

% for continuous predictors
convkernel = nan( length(t0)*length(tau2), length(x));
ctr=0;
for jj=1:length(t0)
for kk=1:length(tau2)
    ctr=ctr+1;
    hh = (x-t0(jj)).*((x-t0(jj))>0);
    convkernel(ctr, :) =   (hh)/tau1 .* exp(-(hh)/tau2(kk));
    convkernel(ctr,:)  =    convkernel(ctr,:)/sum(convkernel(ctr,:));   %normalise 
end
end


%% for puff
dt = 1000/300;
x2 = 0:dt:1000;
t02 = [0 100 200 300];
tau22 = [40 100 200];
convkernel2 = nan( length(t02)*length(tau22), length(x2));
ctr=0;
for jj=1:length(t02)
for kk=1:length(tau22)
    ctr=ctr+1;
    hh = (x2-t02(jj)).*((x2-t02(jj))>0);
    convkernel2(ctr, :) =   (hh)/tau1 .* exp(-(hh)/tau22(kk));
    convkernel2(ctr,:)  =    convkernel2(ctr,:)/sum(convkernel2(ctr,:));   %normalise 
end
end