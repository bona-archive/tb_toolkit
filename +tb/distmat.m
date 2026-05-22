function D = distmat(x, y)
% TB.DISTMAT  Pairwise Euclidean distance matrix.
%
%   D = tb.distmat(x, y)
%
%   x, y: column or row vectors of site coordinates.
%   D(i,j) = distance between site i and site j.

arguments (Input)
    x (:,1) double {mustBeReal}
    y (:,1) double {mustBeReal}
end

D = sqrt((x(:) - x(:)').^2 + (y(:) - y(:)').^2);

end
