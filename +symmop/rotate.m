function [new_basis] = rotate(basis, theta, rotation_center)

arguments (Input)
    basis (:,2)
    theta (1,1)
    rotation_center (1,2)
end

rotation_matrix = [cos(theta), -sin(theta);
                   sin(theta), cos(theta)];
rotation_matrix_transpose = transpose(rotation_matrix);
new_basis = [basis; (basis-rotation_center)*rotation_matrix_transpose + rotation_center];

tol = 1e-8;
new_basis = uniquetol(new_basis, tol, "ByRows",true);

end