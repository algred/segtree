function flag = inBox(C, L)
% INBOX   For a list of points, decide whether they are in BB
%
% input:
% C: 2 * 4 matrix of BB coordinates, in clockwise order
% L: 2 * N matrix of points' coordinates
% 
% output:
% flag: 1 * N vector, flag(i) > 0 if L(i) is in BB

