function E = engTrackSeg(prm)
% TRACKSEG   Using color, motion and spatial priors to produce an energy
% map for tracking a segment
%
% fields of prm:
%   cp, mp, color_map, motion_map, x0, y0, sp_sz, PixelIdxList
%

% compute spatial prior
[rows cols] = size(prm.color_map);
h = fspecial('gaussian', prm.sp_sz * 2, prm.sp_sz);
spatial_pm = zeros(rows, cols);
xidx = max(1, prm.x0 - sp_sz/2) : min(cols, prm.x0 - sp_sz/2);
yidx = max(1, prm.y0 - sp_sz/2) : min(rows, prm.y0 - sp_sz/2);
spatial_pm(yidx,  xidx) = h(yidx - prm.y0 + prm.sp_sz, xidx - prm.x0 + prm.sp_sz);
% color prior
color_pm = medfilt2(reshape(prm.cp(prm.color_map(:)), rows, cols), [5 5], 'symmetric');
% motion prior
motion_pm = medfilt2(reshape(prm.mp(prm.motion_map(:)), rows, cols), [5 5], 'symmetric');
% final energy map
E = ni(color_pm.*motion_pm.*spatial_pm);

end