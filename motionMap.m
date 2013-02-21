function [ M ] = motionMap( u, v )
%MOTIONFEAT computes motion maps for boundary detection

M(:,:,1) = ni(u);
M(:,:,2) = ni(v);
M(:,:,3) = sqrt(u.*u + v.*v);
M(:,:,4) = ni(u ./ (M(:,:,3) + eps));
M(:,:,5) = ni(v ./ (M(:,:,3) + eps));
M(:,:,3) = ni(M(:,:,3));

end

