function img = drawBox(img, boxs, color)
% DRAWBOX   Draws bounding boxes on an image using given color
%
% boxs: 2n * 4 matrix, every two rows contains the coordinates of vertices
% of a bounding box
% color: a tuple (r, g, b) 

[rows cols chs] = size(img);
N = rows * cols;

for i = 1:2:size(boxs,1)
    c = boxs(i : i+1, :);
    [xs ys] = boxBry(c);
    flg = (xs >= 1) & (xs <= cols) & (ys >= 1) & (ys <= rows);
    xs = xs(flg);
    ys = ys(flg);
    pBox = sub2ind([rows, cols], ys, xs);
    img(pBox) = color(1);
    img(pBox + N) = color(2);
    img(pBox + 2 * N) = color(3);
end

end