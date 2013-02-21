function [cpmap mpmap fmap] = estimatePriors3(trees, imgs, motionMag, params)
% ESTIMATEPRIORS  Estimate the foreground color prior (color, motion) of a video
% and background motion prior, and then estimate the forground maps

cols = size(imgs, 2);
rows = size(imgs, 1);
nfms = size(motionMag, 3);
N = cols * rows;
     
SAMPLE_SIZE = 1000;


nr = 0;
try
bidx = [];
sidx = [];
for i = 1:nfms
    if isempty(trees{i})
        continue;
    end
      
    fidx = [];
    decendt = treeDecendt(trees{i});
    for j = 2:length(trees{i})
        if isempty(decendt{j})
            height = 1;
        else
            height = max([trees{i}(decendt{j}).Depth]);
        end
        nsample = (2^height) * SAMPLE_SIZE;
        idx = ceil(rand(1, nsample) * length(trees{i}(j).PixelIdxList));
        sidx = [sidx (i-1)*N + trees{i}(j).PixelIdxList(idx)];
        fidx = union(fidx, trees{i}(j).PixelIdxList);
        nr = nr + 1;
    end
    bidx = [bidx (i-1)*N + setdiff(1:N, fidx)];
end

% color prior
nImgs = size(imgs, 4);
gmm =  gmdistribution.fit(double([imgs(sidx)' imgs(sidx + nImgs*N)' imgs(sidx + nImgs*N*2)']), 10);

% motion prior
bm = motionMag(bidx);
meanBM = mean(bm);
stdBM = std(bm);

cpmap = zeros(rows, cols, nfms);
mpmap = zeros(rows, cols, nfms);
fmap = zeros(rows, cols, nfms);
for j = 1:nfms
    im = imgs(:,:,:,j);
    im = double([im(1:N)' im(N+1:2*N)' im(2*N+1:3*N)']);
    thisCpmap = pdf(gmm, im);
    thisMpmap = 1 - normpdf(motionMag((j-1)*N + 1:j*N), meanBM, stdBM);
    cpmap(:,:,j) = reshape(thisCpmap, rows, cols);
    mpmap(:,:,j) = reshape(thisMpmap, rows, cols);
    fmap(:,:,j) = medfilt2(cpmap(:,:,j).*mpmap(:,:,j), [5 5], 'symmetric');
end

fmap = fmap / (max(fmap(:))+eps);

catch exception
    getReport(exception)
    keyboard;
end

return;


