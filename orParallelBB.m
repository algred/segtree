function C = orParallelBB(S, or)
% orParallelBB   computes a bounding box for point set S which is parallel to or
%
% S: N * 2 matrix recording [x y] of point's coordinates
% or: given orientation which is the angle with X axis counterclock wise 
%
% C: 2 * 4 matrix contains vertex's coordinates in clockwise order

cosr = cos(or);
sinr = sin(or);
xs = S(:,1);
ys = S(:,2);
x0 = mean(xs);
y0 = mean(ys);

xr = (xs - x0) * cosr - sinr * (ys - y0) + x0;
yr = (xs - x0) * sinr + cosr * (ys - y0) + y0;


minxr = min(xr);
maxxr = max(xr);
minyr = min(yr);
maxyr = max(yr);

CR = [minxr minyr; maxxr minyr; maxxr maxyr; minxr maxyr];

cx = round((CR(:,1) - x0) * cosr + sinr * (CR(:,2) - y0) + x0);
cy = round(-(CR(:,1) - x0) * sinr + cosr * (CR(:,2) - y0) + y0);


C = [cx'; cy'];

end