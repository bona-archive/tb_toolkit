function [Hk] = blochham(k_vec, basis, latvecs, hops, onsite)
% input
  % - k_vec — (1×2) k벡터
  % - basis — (n_orb×2) 서브래티스 위치
  % - latvecs — primitive lattice vectors (a1, a2)
  % - hops — (N×4) 테이블 [i, j, dist, t]
  % - onsite — (n_orb×n_orb) onsite 에너지, diagonal

% output
  % - Hk — (n_orb×n_orb) Hermitian 행렬
arguments (Input)
    k_vec (1,2)
    basis (:,2)
    latvecs (2,2)
    hops (:,4)
    onsite (:,:)
end

assert(size(basis,1) == size(onsite,1), 'basis와 onsite 변수의 sublattice(orbita)의 개수가 다릅니다.')
%% basic setup

%nearest neighbor index, orbital number
n_range = -1:1;
n_orb = size(basis,1);

%onsite energy
Hk = onsite;

% primitive lattice vector a1, a2
a1 = latvecs(1,:);
a2 = latvecs(2,:);

%%
for n1 = n_range
    for n2 = n_range
        % primitive displacement
        R = n1*a1 + n2*a2;
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

end