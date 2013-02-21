function [xs ys] = orientLine(rows, cols, seg)

[ys xs] = ind2sub([rows cols], seg.PixelIdxList);
y1 = min(ys);
y2 = max(ys);
if seg.or == pi/2
    ys = [y1:y2]';
    xs = ones(size(ys))*seg.x0;
elseif seg.or == 0
    xs = [min(xs):max(xs)]';
    ys = ones(size(xs))*seg.y0;
else
    x1 = seg.x0 + abs(y1 - seg.y0) / tan(seg.or);
    x2 = seg.x0 - abs(y2 - seg.y0) / tan(seg.or);
    [xs ys] = bresenham(x1,y1,x2,y2);
    flg = (xs >= 1 & xs <= cols);
    xs = xs(flg);
    ys = ys(flg);
end

end