function [sublattice_ind] = basis2sublat_ind(basis)

arguments (Input)
    basis (:,2)
end

i = 1;
sublattice_ind =[];
for basis_ind = 1: size(basis,1)
    sublattice_ind = [sublattice_ind;i];
    i = i+1;
end

end