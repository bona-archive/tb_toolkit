function [new_basis] = inversion(basis, inversion_center)

arguments (Input)
    basis (:,2)
    inversion_center (1,2)
end

inversion_matrix = [-1, 0;
                    0, -1];
new_basis = [basis; (basis - inversion_center) * inversion_matrix + inversion_center];
tol = 1e-8;
new_basis = uniquetol(new_basis,tol,"ByRows",true);
end