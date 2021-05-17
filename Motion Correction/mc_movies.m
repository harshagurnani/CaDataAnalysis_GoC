%% Show Movement Corr Movie

[nX, nY, nT, nROI] = size(raw);

nrow = ceil(sqrt(nROI));
[mov_raw, mov_mc] = deal(zeros( nX*nrow, nY*nrow, nT ));

% Original [X x Y x T x ROI ]
raw( isnan(raw)) = 0;
% Corrected
mc_raw( isnan( mc_raw) ) = 0;

type = 'uint8'; scale = 2^8;

%rescale each 
for roi=1:nROI

    maxF = prctile( reshape( raw(:,:,:,roi), [nX*nY*nT,1] ), 99 );
    minF = 90;%prctile( reshape( raw(:,:,:,roi), [nX*nY*nT,1] ), 5 );
    
    tmp1 = raw(:,:,:,roi);
    tmp1(tmp1<minF) = 0;
    tmp1 = tmp1*scale/maxF;
    
    tmp2 = mc_raw(:,:,:,roi);
    tmp2(tmp2<minF) = 0;
    tmp2 = tmp2*scale/maxF;
    
    row = floor( (roi-1)/nrow );
    col = mod( roi-1,nrow );
    
    mov_raw( row*nX+1:(row+1)*nX, col*nY+1:(col+1)*nY, :) = tmp1;%raw(:,:,:,roi)*scale/maxF;%
    mov_mc( row*nX+1:(row+1)*nX, col*nY+1:(col+1)*nY, :)  = tmp2;%mc_raw(:,:,:,roi)*scale/maxF;%
end

mov_mc(isnan(mov_mc)) = 0;
mov_raw = uint8(mov_raw);
mov_mc = uint8(mov_mc);

%% write tif
imwrite( mov_raw(:,:,1), 'movie_raw.tif');
imwrite(mov_mc(:,:,1), 'movie_corrected.tif');
for jj=2:nT
imwrite(mov_raw(:,:,jj),'movie_raw.tif','WriteMode','append');
imwrite(mov_mc(:,:,jj),'movie_corrected.tif','WriteMode','append');
end

x = mean(mov_mc,3);
x=uint8(x);
imwrite(x, 'mean_mov.tif');
for jj=2:nT
imwrite(x,'mean_mov.tif','WriteMode','append');
end