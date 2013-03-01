function [ M ] = motionMap( u, v )
%MOTIONFEAT computes motion maps for boundary detection

M(:,:,1) = ni(u);
M(:,:,2) = ni(v);
M(:,:,3) = sqrt(u.*u + v.*v);
M(:,:,4) = min(2, max(0, (u ./ (M(:,:,3) + sqrt(2)) + 1))) / 2;
M(:,:,5) = min(2, max(0, (v ./ (M(:,:,3) + sqrt(2)) + 1))) / 2;
M(:,:,3) = ni(M(:,:,3));

end

