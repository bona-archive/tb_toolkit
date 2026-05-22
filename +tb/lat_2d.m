function [x, y] = lat_2d(Nx, Ny, latvec_x, latvec_y, unitcell)
% TB.LAT2D  Generate 2D real-space lattice sites.
%
%   [x, y] = tb.lat2d(Nx, Ny, latvec_x, latvec_y, unitcell)
%
%   Output x, y are (Nx x Ny x sublat_num) arrays.
%   Flatten with x(:), y(:) for use with tb.distmat / tb.drawbonds.

arguments (Input)
    Nx       (1,1) {mustBeInteger, mustBePositive}
    Ny       (1,1) {mustBeInteger, mustBePositive}
    latvec_x (1,2) double {mustBeReal}
    latvec_y (1,2) double {mustBeReal}
    unitcell (:,2) double {mustBeReal}
end
arguments (Output)
    x
    y
end

sublat_num = size(unitcell, 1);
x = zeros(Nx, Ny, sublat_num);
y = zeros(Nx, Ny, sublat_num);

[X_ind, Y_ind] = ndgrid(0:Nx-1, 0:Ny-1);

BaseX = X_ind .* latvec_x(1) + Y_ind .* latvec_y(1);
BaseY = X_ind .* latvec_x(2) + Y_ind .* latvec_y(2);

for s = 1:sublat_num
    x(:,:,s) = BaseX + unitcell(s,1);
    y(:,:,s) = BaseY + unitcell(s,2);
end

end
