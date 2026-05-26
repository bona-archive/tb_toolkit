function [super_basis, super_sublat, prim_super_sublat_map] = lat_2d2superbasis(basis_x,basis_y)

super_basis = [basis_x(:),basis_y(:)];
[~,~,prim_sublat] = ind2sub(size(basis_x),1:numel(basis_x));
prim_sublat = prim_sublat';

super_sublat = [1:numel(basis_x)]';
prim_super_sublat_map = [prim_sublat,super_sublat];
end