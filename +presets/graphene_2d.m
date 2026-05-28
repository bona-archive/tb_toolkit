function [basis,latvec_x,latvec_y,latvecs, a0] = graphene_2d()

% basis
basis = [sqrt(3)/3,0, 1;
    2*sqrt(3)/3, 0, 2];

a0 = sqrt(3)/3;

% lattice vector
latvec_x = [sqrt(3)/2, 1/2];
latvec_y = [sqrt(3)/2, -1/2];
latvecs = [latvec_x;latvec_y];


end