function [new_basis] = rotate(basis, theta, rotation_center,option)

% tayloring시 좌표 충돌을 피할 것

arguments (Input)
    basis (:,3)
    theta (1,1)
    rotation_center (1,2)
    option.IncludeOriginal (1,1) logical = true
end

xy = basis(:,1:2);
s = basis(:,3);

rotation_matrix = [cos(theta), -sin(theta);
                   sin(theta), cos(theta)];

rotation_matrix_transpose = transpose(rotation_matrix);

if option.IncludeOriginal
    new_xy = [xy; (xy-rotation_center)*rotation_matrix_transpose + rotation_center];
    new_s = [s;s];
else
    new_xy = [(xy-rotation_center)*rotation_matrix_transpose + rotation_center];
    new_s = [s];
end

new_basis = [new_xy, new_s];

tol = 1e-8;
new_basis = uniquetol(new_basis, tol, "ByRows",true);

end