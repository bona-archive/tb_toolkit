function [reci_x, reci_y] = prim2reci(prim_x, prim_y)
% TB.PRIM2RECI  Primitive lattice vectors -> reciprocal lattice vectors.
%
%   [reci_x, reci_y] = tb.prim2reci(prim_x, prim_y)
%
%   Satisfies: prim_i · reci_j = 2π δ_{ij}

arguments (Input)
    prim_x (1,2) double {mustBeReal}
    prim_y (1,2) double {mustBeReal}
end

z = [0 0 1];
a1 = [prim_x, 0];
a2 = [prim_y, 0];

area = abs(cross(a1, a2));
area = area(3);

b1 = 2*pi * cross(z, a2) / area;
b2 = 2*pi * cross(a1, z) / area;

reci_x = b1(1:2);
reci_y = b2(1:2);

end
