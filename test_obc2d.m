%% call graphene preset

[prim_basis,latvec_x,latvec_y,latvecs, a0] = presets.graphene_2d();

%% make supercell coords from primitive basis

Ny = 10;
[basis_x, basis_y] = lat.lat_2d(2, Ny, latvec_x, latvec_y, prim_basis);
super_basis = obc.primbasisxy2superbasis(basis_x,basis_y);

%% supercell connectivity check

figure
hold on
[xvec,yvec]=draw.Draw_coupling_graph(basis_x,basis_y,a0,true);

%% k space

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
kpts = 100;
kk = linspace(-pi,pi,kpts);

%% solve eigen value problem for each k

Delta = 0.3;
onsite = [Delta 0;
          0 Delta];

t = -1;
prim_hops = [1,2,a0,t;
        2,1,a0,t];
super_hops = primhop2superhop();

band = zeros(kpts, size(basis_x,2));
eigenstates = zeros(kpts, size(basis_x,2), size(basis_x,2));

for k_ind = 1:length(kk)
    k_vec = kk(k_ind) * reci_y;
    [Hk] = bloch.blochham_2d(k_vec, super_basis, latvecs, hops, onsite);
    [u,d] = eig(Hk);
    band(k_ind,:) = diag(d);
    eigenstates(k_ind,:,:) = u;
end