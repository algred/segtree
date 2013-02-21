function [xs ys] = boxBry(B)
% BOXBRY     Returns the boundary points of a box B.
% 
% B is 2 * 4, each column is [x;y] of a vertex

[xs1 ys1] = bresenham(B(1,1),B(2,1),B(1,2),B(2,2));
[xs2 ys2] = bresenham(B(1,2),B(2,2),B(1,3),B(2,3));
[xs3 ys3] = bresenham(B(1,3),B(2,3),B(1,4),B(2,4));
[xs4 ys4] = bresenham(B(1,4),B(2,4),B(1,1),B(2,1));

xs = [xs1; xs2; xs3; xs4];
ys = [ys1; ys2; ys3; ys4];
end