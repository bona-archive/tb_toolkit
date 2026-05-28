function [new_basis] = mirror(basis, theta, point,option)

% tayloring시 좌표 충돌을 피할 것

arguments (Input)
    basis (:,3)
    theta (1,1)
    point (1,2)
    option.IncludeOriginal (1,1) logical = true
end

xp = point(1);
yp = point(2);

xy = basis(:,1:2);
s = basis(:,3);

xy_tmp = [xy, ones(1,size(xy,1)).'];
mirror_matrix = [cos(2*theta), sin(2*theta), -xp*cos(2*theta)-yp*sin(2*theta)+xp;
                 sin(2*theta), -cos(2*theta), -xp*sin(2*theta)+yp*cos(2*theta)+yp;
                 0           , 0            , 1                                 ];
if option.IncludeOriginal
    new_xy_tmp = [xy_tmp; xy_tmp*mirror_matrix.'];
    new_xy = new_xy_tmp(:, 1:2);
    new_s = [s;s];
else
    new_xy_tmp = [xy_tmp*mirror_matrix.'];
    new_xy = new_xy_tmp(:, 1:2);
    new_s = [s];
end
new_basis = [new_xy, new_s];

tol = 1e-8;
new_basis = uniquetol(new_basis, tol, 'ByRows',true);

end