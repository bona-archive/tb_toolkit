function [new_basis] = inversion(basis, inversion_center, option)

% tayloring시 좌표 충돌을 피할 것


arguments (Input)
    basis (:,3)
    inversion_center (1,2)
    option.IncludeOriginal (1,1) logical = true
end

inversion_matrix = [-1, 0;
                    0, -1];
xy = basis(:,1:2);
s = basis(:,3);

if option.IncludeOriginal
    new_xy = [xy;(xy - inversion_center) * inversion_matrix + inversion_center];
    new_s = [s;s];
else
    new_xy = [(xy - inversion_center) * inversion_matrix + inversion_center];
    new_s = [s];
end
new_basis = [new_xy,new_s];
tol = 1e-8;
new_basis = uniquetol(new_basis,tol,"ByRows",true);
end