function [K_x, K_y] = reci2BZ(reci_x, reci_y, num_x, num_y, options)
% TB.RECI2BZ  Generate 2D k-grid over the first Brillouin zone.
%
%   [K_x, K_y] = tb.reci2BZ(reci_x, reci_y, num_x, num_y)
%   [K_x, K_y] = tb.reci2BZ(..., open_end=true)
%
%   open_end=false (default): linspace(-0.5, 0.5, num) — includes both edges.
%   open_end=true : (0:num-1)/num - 0.5 — excludes right edge, safe for
%                   periodic sampling (Wilson loop, Chern number).

arguments (Input)
    reci_x  (1,2) double {mustBeReal}
    reci_y  (1,2) double {mustBeReal}
    num_x   (1,1) {mustBeInteger, mustBePositive}
    num_y   (1,1) {mustBeInteger, mustBePositive}
    options.open_end (1,1) logical = false
end

if options.open_end
    fx = (0:num_x-1) / num_x - 0.5;
    fy = (0:num_y-1) / num_y - 0.5;
else
    fx = linspace(-0.5, 0.5, num_x);
    fy = linspace(-0.5, 0.5, num_y);
end

[F1, F2] = meshgrid(fx, fy);

K_x = F1 .* reci_x(1) + F2 .* reci_y(1);
K_y = F1 .* reci_x(2) + F2 .* reci_y(2);

end
