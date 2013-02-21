function [cp cpmap mpmap fmap] = estimatePriors2(trees, imgs, motionMag, params)
% ESTIMATEPRIORS  Estimate the foreground color prior (color, motion) of a video
% and background motion prior, and then estimate the forground maps

cols = size(imgs, 2);
rows = size(imgs, 1);
nfms = size(motionMag, 3);
N = cols * rows;
     
CH = zeros(1, params.pr_cnbins);
MH = zeros(1, params.pr_maxbgm);

indC = zeros(rows, cols, nfms);
indM = zeros(rows, cols, nfms);
[~, cmap] = rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,1))), params.pr_cnbins);

nr = 0;
try
bidx = [];
for i = 1:nfms
    if isempty(trees{i})
        continue;
    end
    
    map = rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,i))), cmap) + 1;
    map = medfilt2(map, [3 3], 'symmetric');
    indC(:,:,i) = map;    
    fidx = [];
    decendt = treeDecendt(trees{i});
    for j = 2:length(trees{i})
        if isempty(decendt{j})
            height = 1;
        else
            height = max([trees{i}(decendt{j}).Depth]);
        end
        counts = histc(map(trees{i}(j).PixelIdxList), 1:params.pr_cnbins);
        CH = CH + (counts / max(counts)) * (2^height);
        fidx = union(fidx, trees{i}(j).PixelIdxList);
        nr = nr + 1;
    end
    bidx = [bidx (i-1)*N + setdiff(1:N, fidx)];
end

cp = CH / sum(CH);
bm = motionMag(bidx);
meanBM = mean(bm);
stdBM = std(bm);

cpmap = zeros(rows, cols, nfms);
mpmap = zeros(rows, cols, nfms);
fmap = zeros(rows, cols, nfms);
for j = 1:nfms
    thisIndC = indC(:,:,j);
    thisCpmap = cp(thisIndC(:));
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


