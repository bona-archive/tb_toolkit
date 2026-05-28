clear all
close all
%% build basis with grain boundary

[basis,latvec_x_preset,~,~, a0]=presets.graphene_2d();

basis_1 = basis(1,:);
basis_2 = basis(2,:);
basis_3 = symmop.translate(basis_1,latvec_x_preset,"IncludeOriginal",false);
basis_4 = symmop.translate(basis_2,latvec_x_preset,"IncludeOriginal",false);

basis_5 = basis_1;
basis_5(3) = 2;
basis_6 = basis_2;
basis_6(3) = 1;
basis_7 = basis_3;
basis_7(3) = 2;
basis_8 = basis_4;
basis_8(3) = 1;

basis_A = symmop.merge(basis_1,basis_2,basis_3,basis_4);
basis_B = symmop.merge(basis_5,basis_6,basis_7,basis_8);

N_iter = 10;
basis_B = symmop.translate(basis_B, (N_iter+1)*[sqrt(3),0],"IncludeOriginal",false);
%% build supercell with grain boundary

latvec_x = [sqrt(3), 0];
latvec_y = [0, -1];
latvecs = [latvec_x;latvec_y];

basis_grainA = basis_A;
for i = 1: N_iter
    basis_grainA = symmop.translate(basis_grainA,latvec_x, "IncludeOriginal", true); 
end

basis_grainB = basis_B;
for i = 1: N_iter
    basis_grainB = symmop.translate(basis_grainB, latvec_x, "IncludeOriginal", true);
end

basis_AB = symmop.merge(basis_grainA, basis_grainB);
%% test draw supercell

figure
draw.Draw_coupling_graph(basis_AB(:,1), basis_AB(:,2), a0,'k', true);

%% test draw bulk (latvec_pbc direction)

figure
test_iter = 30;
test_basis = basis_AB;
latvec_pbc = latvec_y;
for iter = 1: test_iter
    test_basis = symmop.translate(test_basis,latvec_pbc);
end
draw.Draw_coupling_graph(test_basis(:,1),test_basis(:,2),a0,'k',true);
%% build hop, onsite for supercell

Delta1 = 0;
Delta2 = 2;
Delta_vec = [Delta1, Delta2];
onsite = bloch.build_onsite(basis_AB, Delta_vec);

t = -1;
hops = [1,2,a0,t;
        2,1,a0,t];

%% build bloch Hamiltonian

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
k_pts = 100;
kk = linspace(-1/2,1/2,k_pts);

band = zeros(k_pts, size(basis_AB,1));
for k_ind = 1: length(kk)
    k_vec = kk(k_ind)*reci_y;
    [Hk] = bloch.blochham_2d(k_vec,basis_AB,latvec_pbc, hops, onsite);
    [u,d] = eig(Hk);
    band(k_ind,:)=diag(d);
end

figure
hold on
plot(kk, band)