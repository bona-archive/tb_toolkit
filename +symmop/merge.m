function [merged_basis] = merge(varargin)

arguments (Repeating)
    varargin
end

merged_basis = vertcat(varargin{:});
tol = 1e-8;
merged_basis = uniquetol(merged_basis, tol, 'ByRows',true);
end