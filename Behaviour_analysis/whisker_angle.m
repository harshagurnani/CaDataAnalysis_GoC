% from tracking

%% set whisker pad vector ( line with which angle is measured )


%% loaded whiskers - w1, w2, w3
% allw = dir( [cpath, '*filtered*.csv']);

%% get angle
nTm = length( whiskers{1} );
nWhisk = length( whiskers );

vm = nan( nTm, 2, nWhisk );
angle = nan( nTm, nWhisk );
s2 = [0, wpad_vm(1); 0, wpad_vm(2) ];

for ww = 1:nWhisk
for jj = 1:nTm
    vm(jj,:,ww) = polyfit(whiskers{ww}(jj,[1,3,5]),whiskers{ww}(jj,[2,4,6]),1);
    s1 = [0, 1; 0, vm(jj,1,ww)];

    angle(jj,ww) = subspace( s1, s2 );
end
end