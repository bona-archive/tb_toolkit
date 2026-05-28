function [new_basis] = translate(basis, translation_vector,option)

% tayloring시 좌표 충돌을 피할 것

arguments (Input)
    basis (:,3)
    translation_vector (1,2)
    option.IncludeOriginal (1,1) logical = true
end

xy = basis(:,1:2);
s = basis(:,3);

if option.IncludeOriginal
    new_xy = [xy; xy+translation_vector];
    new_s = [s;s];
else
    new_xy = [xy+translation_vector];
    new_s = [s];
end

new_basis = [new_xy,new_s];

tol = 1e-8;
new_basis = uniquetol(new_basis, tol, 'ByRows', true);
end