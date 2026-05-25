function [Hk] = blochham_2d(k_vec, basis, latvecs, hops, onsite)
% 2차원 공간에서의 블로흐 해밀토니안 생성. 2D 벌크 블로흐는 latvecs에 basis 둘 다 넣을 것
% OBC를 주고 싶은 basis vector는 latvecs에 넣지 않으면 됨. 
% input
  % - k_vec — (1×2) k벡터
  % - basis — (n_orb×2) 서브래티스 위치
  % - latvecs — PBC 조건을 넣을 primitive lattice vectors (a1, a2)로 주면 2D bulk,
  % 하나만 주면 그 방향으로만 PBC.
  % - hops — (N × 4) 테이블 [i, j, dist, t]
  % - onsite — (n_orb × n_orb) onsite 에너지, diagonal
% output
  % - Hk — (n_orb×n_orb) Hermitian 행렬
arguments (Input)
    k_vec (1,2)
    basis (:,2)
    latvecs (:,2)
    hops (:,4)
    onsite (:,:)
end

assert(size(basis,1) == size(onsite,1), 'basis와 onsite 변수의 sublattice(orbita)의 개수가 다릅니다.')
%% basic setup
% PBC를 줄 차원 개수
D = size(latvecs,1);

% D가 결정되면 Nearest Neighbor cell로 3^D사이즈의 grid 동적 생성
n_range = -1:1;

args = repmat({n_range},1,D);
grids = cell(1,D);
[grids{:}] = ndgrid(args{:});
result = cellfun(@(g) g(:), grids, 'UniformOutput',false);
NNcell_grid = [result{:}];

% sublattice(orbit) 개수
n_orb = size(basis,1);

%onsite energy
Hk = onsite;

% primitive lattice vector a1, a2
% a1 = latvecs(1,:);
% a2 = latvecs(2,:);

%%
for grid_ind = 1: size(NNcell_grid, 1)
    grid_coords = NNcell_grid(grid_ind,:);
    R = zeros(1, 2);
    for coord_ind = 1: size(NNcell_grid,2)
        % primitive displacement
        R = R + grid_coords(coord_ind)*latvecs(coord_ind,:);
    end
        for i = 1: n_orb
            for j = 1: n_orb
                % sublattice displacement
                b1 = basis(i,:);
                b2 = basis(j,:);
                d_vec = R + b2 - b1;
                dist = norm(d_vec);

                % encoding bloch phase
                match = hops(:,1) == i & hops(:,2) == j & abs(hops(:,3)-dist)<1e-6;
                if any(match)
                    Hk(i,j) = Hk(i,j) + hops(match, 4) * exp(1i * dot(k_vec,d_vec));
                end
            end
        end

end

end