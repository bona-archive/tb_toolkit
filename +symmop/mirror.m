function [new_basis] = mirror(basis, theta, point)

arguments (Input)
    basis (:,2)
    theta (1,1)
    point (1,2)
end

xp = point(1);
yp = point(2);

basis_tmp = [basis, ones(1,size(basis,1)).'];
mirror_matrix = [cos(2*theta), sin(2*theta), -xp*cos(2*theta)-yp*sin(2*theta)+xp;
                 sin(2*theta), -cos(2*theta), -xp*sin(2*theta)+yp*cos(2*theta)+yp;
                 0           , 0            , 1                                 ];
new_basis_tmp = [basis_tmp; basis_tmp*mirror_matrix.'];
tol = 1e-8;
new_basis_tmp = uniquetol(new_basis_tmp, tol, 'ByRows',true);
new_basis = new_basis_tmp(:, 1:2);
end