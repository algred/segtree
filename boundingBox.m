function C = boundingBox(xs, ys)
% BOUNDINGBOX  Compute a bounding box for a set of points. The box is
% aligned to be parallel with the axis of least inertia
%
% xs, ys: 1 * N vector contains the coordinates of points
% C: 2 * 4 matrix where each column is coordinate [x; y] of a bounding box
% vertex and the vertices are listed in clockwise order

[x0 y0 or] = segProperty([xs' ys']);
C = orParallelBB([xs' ys'], or);

end