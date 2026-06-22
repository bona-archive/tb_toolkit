clear all
close all
%% 

load("+custom_colorbar/CustomColormap.mat")
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

%% build hop, onsite for supercell

Delta1 = 5;
Delta2 = 0;
Delta_vec = [Delta1, Delta2];
onsite = bloch.build_onsite(basis_AB, Delta_vec);

t = 0;
t1 = 0.5;
a1 = norm(latvec_x_preset);
hops = [1,2,a0,t;
        2,1,a0,t;
        2,2,a0,t;
        1,1,a1,t1;
        2,2,a1,t1];
% hops = [1,2,a0,t;
%         2,1,a0,t;
%         2,2,a0,t];
% hops = [
%         1,1,a1,t1;
%         2,2,a1,t1];
%% test draw supercell

figure
draw.Draw_coupling_graph(basis_AB(:,1), basis_AB(:,2), basis_AB(:,3), hops, true,basis_AB(:,3));

%% find grain boundary index using mask

epsilon = 0.2;

left_bot = [max(basis_grainA(:,1))-epsilon min(basis_grainA(:,2))-epsilon];

left_top = [max(basis_grainA(:,1))-epsilon max(basis_grainA(:,2))+epsilon];

right_bot = [min(basis_grainB(:,1))+epsilon min(basis_grainB(:,2))-epsilon];

right_top = [min(basis_grainB(:,1))+epsilon max(basis_grainB(:,2))+epsilon];

corners = [left_bot;
           left_top;
           right_top;
           right_bot;
           left_bot];
figure
hold on
draw.Draw_coupling_graph(basis_AB(:,1), basis_AB(:,2), basis_AB(:,3), hops, true,basis_AB(:,3));
plot(corners(:,1), corners(:,2))

in = inpolygon(basis_AB(:,1),basis_AB(:,2),corners(1:4,1),corners(1:4,2));
in_x = basis_AB(in,1);
in_y = basis_AB(in,2);
gb_idx = find(in == 1);
%% test draw bulk (latvec_pbc direction)

figure
test_iter = 30;
test_basis = basis_AB;
latvec_pbc = latvec_y;
for iter = 1: test_iter
    test_basis = symmop.translate(test_basis,latvec_pbc);
end
draw.Draw_coupling_graph(test_basis(:,1),test_basis(:,2),test_basis(:,3), hops,true,test_basis(:,3));
% colormap([1,0,0;0,0,1])
% draw.Draw_coupling_graph(test_basis(:,1),test_basis(:,2),1,'k',true,test_basis(:,3));

%% build bloch Hamiltonian

[reci_x,reci_y]=to.prim2reci(latvec_x,latvec_y);
k_pts = 100;
kk = linspace(-1/2,1/2,k_pts);

band = zeros(k_pts, size(basis_AB,1));
eigenstate = zeros(k_pts, size(basis_AB,1), size(basis_AB,1));
Wgb = zeros(k_pts,size(basis_AB,1));

for k_ind = 1: length(kk)
    k_vec = kk(k_ind)*reci_y;
    [Hk] = bloch.blochham_2d(k_vec,basis_AB,latvec_pbc, hops, onsite);
    [u,d] = eig(Hk);
    band(k_ind,:)=diag(d);
    eigenstate(k_ind,:,:) = u;
    Wgb(k_ind,:) = sum(abs(u(gb_idx,:)).^2, 1);
end

figure
hold on
plot(kk, band)

%%

% grain boundary amplitude를 색깔로 나타낸 full band 플롯

figure
hold on
X = repmat(kk',1,size(band,2));
X = X(:);
scatter(X, band(:), 12, Wgb(:), 'filled');
colormap(CustomColormap)

colorbar; xlabel('k_y'); ylabel('E_{k_y}');
%% 

% grain boundary amplitude로 선별된 밴드
edge_inds = find(sum(Wgb,1)>mean(sum(Wgb,1))+3*std(sum(Wgb,1)));

for edge_ind = edge_inds
    edge_band1 = band(:,edge_ind);
    % edge_band2 = band(:,end);
    edge_band_diff = diff(edge_band1);
    tol_diff = 0.0001;
    diff0_ind = find(abs(edge_band_diff) < tol_diff)
    
    % scanning omega 선별
    omega0 = edge_band1(diff0_ind);
    
    % LDOS 계산
    eta = 1e-5;
    rho = zeros(size(basis_AB,1),1);
    for k_ind = 1:k_pts
        u_tmp = squeeze(eigenstate(k_ind,:,:));
        d_tmp = squeeze(band(k_ind,:));
        rho = rho + sum(abs(u_tmp).^2 .* (-1/pi)*imag(1./(omega0-diag(d_tmp).'+1i*eta)), 2);
    end
    
    rho = rho/k_pts;
    
    
    figure
    hold on
    % draw.Draw_coupling_graph(basis_AB(:,1), basis_AB(:,2), basis_AB(:,3),hops,true,basis_AB(:,3));
    
    scatter(basis_AB(:,1), basis_AB(:,2), 40, rho, 'filled')
    colorbar
    colormap(CustomColormap)
    axis equal off
    plot(basis_AB(:,1),rho/max(abs(rho))+1,'k')
end
