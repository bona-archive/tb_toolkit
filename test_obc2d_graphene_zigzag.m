clear all
% close all

%% call graphene preset

[prim_basis,latvec_x,latvec_y,latvecs, a0] = presets.graphene_2d();

%% make supercell coords from primitive basis
Nx = 1;
Ny = 30;
latvec_pbc = latvec_x;
[superlat_2d_x, superlat_2d_y] = lat.lat_2d(Nx, Ny, latvec_x, latvec_y, prim_basis);
[super_basis, super_sublat, prim_super_sublat_map] = obc.lat_2d2superbasis(superlat_2d_x,superlat_2d_y);

superlatvec_x = Nx*latvec_x;
superlatvec_y = Ny*latvec_y;
super_latvecs = [superlatvec_x;superlatvec_y];
% [super_basis_x, super_basis_y] = lat.lat_2d(1, 1, Nx*latvec_x, Ny*latvec_y, super_basis);

%% supercell connectivity check

figure
hold on
[xvec,yvec]=draw.Draw_coupling_graph(superlat_2d_x,superlat_2d_y,a0,'k',true);

%% k space

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
kpts = 100;
kk = linspace(-0.5,0.5,kpts);

%% solve eigen value problem for each k

Delta = 0;
onsite_tmp = Delta*diag(ones(1,Nx*Ny*2));

t = -1;
prim_hops = [1,2,a0,t;
        2,1,a0,t];

% [super_hops] = obc.primhop2superhop(prim_hops,superlat_2d_x,superlat_2d_y, latvec_pbc);

band = zeros(kpts, size(superlat_2d_x,2)+size(superlat_2d_y,2));
eigenstates = zeros(kpts, size(superlat_2d_x,2)+size(superlat_2d_y,2), size(superlat_2d_x,2)+size(superlat_2d_y,2));

for k_ind = 1:length(kk)
    k_vec = kk(k_ind) * reci_x;
    [Hk] = bloch.blochham_2d(k_vec, super_basis, latvec_pbc, prim_hops, onsite_tmp);
    [u,d] = eig((Hk+Hk')/2);
    band(k_ind,:) = diag(d);
    eigenstates(k_ind,:,:) = u;
end

% Plot the band structure
figure;
plot(kk, band);
xlabel('k-point');
ylabel('Energy');
title('Band Structure of Graphene Supercell');
grid on;