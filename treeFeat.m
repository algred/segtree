function tree = treeFeat(tree, img, params)
[rows cols chs] = size(img);
N = rows * cols;
gray = rgb2gray(img);

for i = 1:length(tree)
    sz = length(tree(i).PixelIdxList);
    m = zeros(5 + params.bw_quant_num, sz);
    
    % color
    m(1,:) = img(tree(i).PixelIdxList);
    m(2,:) = img(tree(i).PixelIdxList + N);
    m(3,:) = img(tree(i).PixelIdxList + 2 * N);
    
    % gradient
    [gx gy] = gradient(double(gray));
    m(4,:) = gx(tree(i).PixelIdxList);
    m(5,:) = gy(tree(i).PixelIdxList);
    
    % intensity features
    f = getFrameExpBlur(gray, params.bw_quant_num, params.bw_sigma_blur);
    for j = 1:params.bw_quant_num
        m(5 + j, :) = f(tree(i).PixelIdxList + (j-1) * N);
    end
    
    tree(i).Feat = m;
end


end