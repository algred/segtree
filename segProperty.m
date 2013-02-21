function [x y or] = segProperty(S)
%SEGPROPERTY  computes the center and orientation of segment S
%
% S is a N * 2 matrix containing the coordinates of points in the segment

S(:,2) = -S(:,2);

x = mean(S(:,1));
y = mean(S(:,2));


a = sum((S(:, 1) - x) .* (S(:, 1) - x));
b = 2 * sum((S(:, 1) - x) .* (S(:, 2) - y));
c = sum((S(:, 2) - y) .* (S(:, 2) - y));

if b == 0 && a == c
    or = 0;
else
    d = a - c;
    sin2 = b / sqrt(b * b + d * d);
    cos2 = d / sqrt(b * b + d * d);
    if sin2 == 0
    	or = pi / 4;
    elseif cos2 == 0
    	or = pi / 2;
    elseif sin2 > 0 && cos2 > 0
    	or = asin(sin2) / 2;
    elseif sin2 > 0 && cos2 < 0
    	or = acos(cos2) / 2;
    elseif sin2 < 0 && cos2 < 0
    	or = acos(-cos2)/2+pi/2;
    elseif sin2 < 0 && cos2 > 0
    	or = pi - asin(-sin2)/2;
    end
end

y = -y;

end