function [new_basis] = translate(basis, translation_vector)

arguments (Input)
    basis (:,2)
    translation_vector (1,2)
end

new_basis = [basis; basis+translation_vector];

tol = 1e-8;
new_basis = uniquetol(new_basis, tol, 'ByRows', true);
end