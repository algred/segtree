% author: Marius Leordeanu
% last modified: October 2012

%--------------------------------------------------------------------------
% Gb Code Version 1
%--------------------------------------------------------------------------

% INPUT:  I - rgb color image
%         M - motion maps: u, v, normalized u, normalized v, magnitude
%         P - foreground prior map

%--------------------------------------------------------------------------

% OUTPUT: 

%         gb_thin_CSG  - thin Gb boundaries using Color (C), Soft-segmentation (S) and Geometric grouping (G) 
%                        with nonlocal maxima suppression 
%         gb_thin_CS   - thin Gb boundaries using Color (C) and Soft-segmentation (S) 
%                        with nonlocal maxima suppression 
%         gb_CS        - continuous Gb boundaries using Color (C) and Soft-segmentation (S) 
%                        without nonlocal maxima suppression 
%
%         orC          - orientation computed from color channels, in
%                        degrees between 0 and pi
%
%         edgeImage      - color image with different pieces of contours in
%                          different colors; use it for visualization
%
%         edgeComponents - of the same size as I, with each contour having
%                          a unique integer assigned; a possible use is to find all
%                          pixels belonging to one contour 


% NOTE: gb_thin_CSG has better accuracy than gb_thin_CS

% This code is for research use only. 
% It is based on the following paper, which should be cited:

%  Marius Leordeanu, Rahul Sukthankar and Cristian Sminchisescu, 
% "Efficient Closed-form Solution to Generalized Boundary Detection", 
%  in ECCV 2012
%
%
%  Modified by Shugao Ma

function [gb_CSMPG, orC] = Gb_CSMPG(I, M, P)

[nRows, nCols, aux] = size(I);
imDiag = norm([nRows, nCols]);

wS = round(0.041*imDiag);
wC = round(0.016*imDiag);
wM = wC;

alpha_AB = 1.9;
alpha_M = 0.5;
alpha_P = 2;

disp('soft-segmentation ...');

seg = softSegs(I);

disp('Gb: color + motion + foreground prior + soft-segmentation ... ');


I = double(rgb2lab(I));
sz = size(I, 3);


for i = 1:size(I, 3)
    D(:,:,i) = ni(I(:,:,i));
end
D(:,:,2:3) = D(:,:,2:3) * alpha_AB;

for i = 1:size(M, 3)
    D(:,:,i+sz) = M(:,:,i) * alpha_M;
end

sz = size(D, 3);
D(:,:,sz+1) = P * alpha_P;

[gbS, orS] =  Gb_data_lambda(seg, wS);
[gbD, orC] = Gb_data_lambda(D, wM);
gb_CSMPG = ni(gbD.*ni(gbS));
orC = orC * pi / 180;

end