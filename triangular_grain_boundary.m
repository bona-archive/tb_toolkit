clear all
close all
%% build super cell with grain boundary

[basis,latvec_x,latvec_y,latvecs, a0]=presets.triangular_2d();

basis = [basis; [1/2,sqrt(3)/2, 2]];
basis_grainA = basis;
iteration_grainA = 10;
for i = 1: iteration_grainA
    basis_grainA = symmop.translate(basis_grainA,[0,sqrt(3)]); 
end

grain_seperation = 0.5;
translate_vector_grainA2B = iteration_grainA*[0,sqrt(3)]+[0,sqrt(3)/2]+[0,grain_seperation];
basis_grainAB = symmop.translate(basis_grainA, translate_vector_grainA2B);
dist_ABcoupling = norm(basis_grainAB(iteration_grainA+2,1:2)-basis_grainAB(3*(iteration_grainA+1),1:2));
figure
draw.Draw_coupling_graph(basis_grainAB(:,1), basis_grainAB(:,2),a0,'k', true);
draw.Draw_coupling_graph(basis_grainAB(:,1), basis_grainAB(:,2),dist_ABcoupling,'r', true);

% test draw
figure
test_iter = 30;
test_basis = basis_grainAB;
for iter = 1: test_iter
    test_basis = symmop.translate(test_basis,latvec_x);
end
draw.Draw_coupling_graph(test_basis(:,1),test_basis(:,2),a0,'k',true);
draw.Draw_coupling_graph(test_basis(:,1),test_basis(:,2),dist_ABcoupling,'r',true);
%% build hop, onsite for supercell

Delta1 = 0;
Delta2 = 0;
Delta_vec = [Delta1, Delta2];
onsite = bloch.build_onsite(basis_grainAB,Delta_vec);

t = -0.5;
t_AB = -0.5;
hops = [1,2,a0,t;
        2,1,a0,t;
        1,2,dist_ABcoupling,t_AB;
        2,1,dist_ABcoupling,t_AB];

%% build bloch Hamiltonian

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
k_pts = 100;
kk = linspace(-1/2,1/2,k_pts);

latvec_pbc = latvec_x;

band = zeros(k_pts, size(basis_grainAB,1));
for k_ind = 1: length(kk)
    k_vec = kk(k_ind)*reci_x;
    [Hk] = bloch.blochham_2d(k_vec,basis_grainAB,latvec_pbc, hops, onsite);
    [u,d] = eig(Hk);
    band(k_ind,:)=diag(d);
end

figure
hold on
plot(kk, band)