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
t_AB = -0.7;
hops = [1,2,a0,t;
        2,1,a0,t;
        1,2,dist_ABcoupling,t_AB;
        2,1,dist_ABcoupling,t_AB];

%% grain index 추출

distmat = lat.distmat(basis_grainAB(:,1), basis_grainAB(:,2));
tol_dist = 1e-8;
I=find(abs(distmat-dist_ABcoupling)<tol_dist);

grain_ind_list = [];
for I_ind = 1: length(I)
    [row,col] = ind2sub(size(distmat),I(I_ind));
    grain_ind_list = [grain_ind_list;row;col];
end

grain_ind_list = unique(grain_ind_list);

gb_idx = grain_ind_list;
%% build bloch Hamiltonian

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
k_pts = 100;
kk = linspace(-1/2,1/2,k_pts);

latvec_pbc = latvec_x;
eigenstate = zeros(k_pts, size(basis_grainAB,1), size(basis_grainAB,1));
band = zeros(k_pts, size(basis_grainAB,1));
Wgb = zeros(k_pts,size(basis_grainAB,1));
for k_ind = 1: length(kk)
    k_vec = kk(k_ind)*reci_x;
    [Hk] = bloch.blochham_2d(k_vec,basis_grainAB,latvec_pbc, hops, onsite);
    [u,d] = eig(Hk);
    eigenstate(k_ind,:,:) = u;
    band(k_ind,:)=diag(d);
    Wgb(k_ind,:) = sum(abs(u(gb_idx,:)).^2, 1);
end

% full band 플롯
figure
hold on
plot(kk, band)

% grain boundary amplitude를 색깔로 나타낸 full band 플롯
figure
hold on
X = repmat(kk',1,size(band,2));
X = X(:);
scatter(X, band(:), 12, Wgb(:), 'filled');
colorbar; xlabel('k'); ylabel('E');

% grain boundary amplitude로 선별된 밴드
edge_band1 = band(:,1);
edge_band2 = band(:,end);
edge_band_diff = diff(edge_band1);
tol_diff = 0.001;
diff0_ind = find(abs(edge_band_diff) < tol_diff);

% scanning omega 선별
omega0 = edge_band1(diff0_ind);

% LDOS 계산
eta = 1e-5;
rho = zeros(size(basis_grainAB,1),1);
for k_ind = 1:k_pts
    u_tmp = squeeze(eigenstate(k_ind,:,:));
    d_tmp = squeeze(band(k_ind,:));
    rho = rho + sum(abs(u_tmp).^2 .* (-1/pi)*imag(1./(omega0-diag(d_tmp).'+1i*eta)), 2);
end

rho = rho/k_pts;


figure
hold on
[xvec,yvec] = draw.Draw_coupling_graph(basis_grainAB(:,1), basis_grainAB(:,2), a0, 'k',false);
line(xvec,yvec,'Color','k')
[xvec,yvec] = draw.Draw_coupling_graph(basis_grainAB(:,1), basis_grainAB(:,2), dist_ABcoupling,'r',false);
line(xvec,yvec,'Color','r')

scatter(basis_grainAB(:,1), basis_grainAB(:,2), 40, rho, 'filled')
colorbar
axis equal
