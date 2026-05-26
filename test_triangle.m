clear all
% close all

%% triangular lattice bulk

% call presets
[basis,latvec_x,latvec_y,latvecs, a0] = presets.triangular_2d();

% make supercell coords from primitive basis
Nx = 10;
Ny = 10;
[x, y] = lat.lat_2d(Nx, Ny, latvec_x, latvec_y, basis);
figure
[xvec,yvec] = draw.Draw_coupling_graph(x,y,a0,true);

% onsite potential
onsite = diag([0,0]);

% hopping dist
neighbor_order = [1];
t12 = -1;
t23 = -1;
t31 = -1;
hops = [1, 2, a0, t12;
        2, 3, a0, t23;
        3, 1, a0, t31]; % [i, j, dist, t]

% high symm points
kpts = {[0,0],[1/sqrt(3), 1/3],[1/2/sqrt(3),1/2],[0,0]};
labels = {'\Gamma','K','M','\Gamma'};

% making k space
N_pts = 100;
[k_path, k_dist, ticks, labels] = bloch.kpath(kpts, labels, N_pts);

% solve eigenvalue
band = zeros(size(basis,1),size(k_path,1));
for k_ind = 1:size(k_path,1)
    k_vec = k_path(k_ind,:);
    Hk = bloch.blochham_2d(k_vec, basis, latvecs, hops, onsite);
    [u,d] = eig((Hk+Hk')/2);
    band(:,k_ind) = sort(diag(d));
end

figure
hold on
for band_ind = 1:size(band,1)
    plot(k_dist,band(band_ind,:));

end
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'FontSize', 30);