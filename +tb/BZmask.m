function [in_BZ, vx, vy] = BZmask(Kx, Ky, vx, vy)
% TB.BZMASK  Mask k-points inside a polygonal Brillouin zone.
%
%   [in_BZ, vx, vy] = tb.BZmask(Kx, Ky, vx, vy)
%
%   vx, vy: polygon vertex coordinates (e.g. from tb.hexBZverts).
%   in_BZ: logical array, same size as Kx.
%
%   Usage after surf plot:
%     [in_BZ, vx, vy] = tb.BZmask(Kx, Ky, vx, vy);
%     EnergyData(~in_BZ) = NaN;

arguments (Input)
    Kx (:,:) double {mustBeReal}
    Ky (:,:) double {mustBeReal}
    vx (1,:) double {mustBeReal}
    vy (1,:) double {mustBeReal}
end

in_BZ = inpolygon(Kx, Ky, vx, vy);

end
