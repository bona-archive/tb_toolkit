function [vx, vy] = hexBZverts(reci_x, reci_y)
% TB.HEXBZVERTS  Compute hexagonal BZ vertices from reciprocal lattice vectors.
%
%   [vx, vy] = tb.hexBZverts(reci_x, reci_y)
%
%   Returns 6 vertex coordinates of the 1st BZ hexagon.
%   Pass directly to tb.BZmask.

arguments (Input)
    reci_x (1,2) double {mustBeReal}
    reci_y (1,2) double {mustBeReal}
end

% BZ corners = (±b1 ± b2)/3 combinations for hexagonal lattice
corners = [
     reci_x + reci_y;
     reci_x;
     reci_y;
    -reci_x - reci_y;
    -reci_x;
    -reci_y
] / 3;  % Wigner-Seitz convention for hex lattice: corners at (b1+b2)/3 etc.

% sort by angle to get convex polygon
angles = atan2(corners(:,2), corners(:,1));
[~, idx] = sort(angles);
corners = corners(idx, :);

vx = corners(:,1)';
vy = corners(:,2)';

end
