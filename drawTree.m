function img = drawTree(img, tree)

try
[rows cols chs] = size(img);
N = rows * cols;

for i = 2:length(tree)
    c = tree(i).BoundingBox;
    [xs ys] = boxBry(c);
    flg = (xs >= 1) & (xs <= cols) & (ys >= 1) & (ys <= rows);
    xs = xs(flg);
    ys = ys(flg);
    pBox = sub2ind([rows, cols], ys, xs);
    img(pBox) = 255;
    img(pBox + N) = 0;
    img(pBox + 2 * N) = 0;
end

catch exception
    getReport(exception)
    keyboard;
end
   

end