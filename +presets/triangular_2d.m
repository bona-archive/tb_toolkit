function [basis,latvec_x,latvec_y,latvecs, a0] = triangular_2d()

% basis
basis = [0, 0];

a0 = 1;

% lattice vector
latvec_x = [1, 0];
latvec_y = [-1/2, sqrt(3)/2];
latvecs = [latvec_x;latvec_y];

end