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
[xvec,yvec] = draw.Draw_coupling_graph(x,y,a0,'k',true);

% onsite potential
onsite = diag([0]);

% hopping dist
neighbor_order = [1];
t = -1;
hops = [1, 1, a0, t]; % [i, j, dist, t]

% high symm points
kpts = {[0,0],[1/2/sqrt(3),1/2],[1/sqrt(3), 1/3],[0,0]};
labels = {'\Gamma','M','K','\Gamma'};

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

%%

kx = linspace(-1/2,1/2,N_pts);
ky = kx;
band = zeros(size(basis,1),N_pts,N_pts);
for kx_ind = 1:size(ky,2)
    for ky_ind = 1:size(ky,2)
        k_vec = [kx(kx_ind),ky(ky_ind)];
        Hk = bloch.blochham_2d(k_vec, basis, latvecs, hops, onsite);
        [u,d] = eig((Hk+Hk')/2);
        band(:,kx_ind,ky_ind) = sort(diag(d));
    end
end
figure
hold on
for band_ind = 1:size(band,1)
    surf(kx,ky,squeeze(band(band_ind,:,:)))
end
view(3)
figure
[kx_mesh,ky_mesh] = meshgrid(kx,ky);
surf(kx_mesh,ky_mesh,2*t*(cos(kx_mesh)+2*cos(1/2*kx_mesh).*cos(sqrt(3)/2*ky_mesh)))