%% graphene

%unit cell
basis = [sqrt(3)/3,0;
    2*sqrt(3)/3, 0];
a0 = 1;

% lattice vector
a1 = [sqrt(3)/2, 1/2];
a2 = [sqrt(3)/2, -1/2];
latvecs = [a1;a2];

% onsite potential
onsite = diag([0,0]);

% hopping dist
neighbor_order = [1];
t = -1;
hops = [1, 2, a0, t]; % [i, j, dist, t]

% high symm points
kpts = {[0,0],[2*pi/sqrt(3), 2*pi/3],[pi/sqrt(3),pi],[0,0]};
labels = {'\Gamma','K','M','\Gamma'};

%% triangle

%unit cell
basis = [[0, 0];              % Site A (Sublattice 1)
[1/2, sqrt(3)/6];    % Site B (Sublattice 2)
[0, sqrt(3)/3];];      % Site C (Sublattice 3)
a0 = 1/sqrt(3);
% lattice vector
a1 = [1, 0];
a2 = [1/2, sqrt(3)/2];
latvecs = [a1;a2];

% onsite potential
Delta = 0.3;
onsite = diag([Delta, 0 , -Delta]);

% hopping dist
neighbor_order = [1];
t_ab = -1;
t_bc = -0.7;
t_ca = -0.2;
hops = [1, 2, a0, t_ab;
        2, 3, a0, t_bc;
        3, 1, a0, t_ca]; % [i, j, dist, t]

% high symm points
kpts = {[0,0],[4*pi/3, 0],[pi, pi/sqrt(3)],[0,0]};
labels = {'\Gamma','K','M','\Gamma'};
%%
% hops = tb.make_hops(basis, latvecs, neighbor_order, t_vals, doPlot=true);


N_pts = 100;
[k_path, k_dist, ticks, labels] = bloch.kpath(kpts,labels, N_pts);

band = zeros(size(basis,1),size(k_path,1));
for k_ind = 1:size(k_path,1)
    k_vec = k_path(k_ind,:);
    Hk = bloch.blochham(k_vec, basis, latvecs, hops, onsite);
    [u,d] = eig((Hk+Hk')/2);
    band(:,k_ind) = sort(diag(d));
end

figure
hold on
for band_ind = 1:size(band,1)
    plot(k_dist,band(band_ind,:));

end
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'FontSize', 30);
