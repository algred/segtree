function a = rotateAngle(C)
% ROTATEBOX    Computes a rotation angle given BB coordinates. 
% This rotation angle will rotate the BB to upright position, 
% i.e the longer edge will be vertical after rotation
%
% Note: the points are listed in clockwise order in C
%       the computed angle is in counter-clockwise

% the lowest point
[~, id1] = max(C(2,:));

% next point in clockwise order
if id1 == 4
    id2 = 1;
else
    id2 = id1 + 1;
end

% previous point in clockwise order
if id1 == 1
    id3 = 4;
else
    id3 = id1 - 1;
end

x1 = C(1, id1) - C(1, id2);
y1 = C(2, id1) - C(2, id2);

if x1 == 0
    v1 = y1 * y1;
else
    v1 = x1 * x1 + y1 * y1;
end

x2 = C(1, id3) - C(1, id1);
y2 = C(2, id1) - C(2, id3);

if y2 == 0
    v2 = x2 * x2;
else
    v2 = x2 * x2 + y2 * y2;
end

if v1 > v2
    if x2 == 0
        a = pi / 2;
    else
        a = -atan( y2 / x2 );
    end
else
    if x1 == 0
        a = pi / 2;
    else
        a = atan( y1 / x1);
    end
end
    
   
    






