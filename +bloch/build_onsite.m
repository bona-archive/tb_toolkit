function [onsite] = build_onsite(basis,Delta_vec)

arguments (Input)
    basis (:,3)
    Delta_vec (1,:)
end
assert(max(basis(:,3)) <= numel(Delta_vec), ...
         'Delta_vec 길이가 sublattice 개수보다 짧습니다.')

Delta_diag = Delta_vec(basis(:,3));
onsite = diag(Delta_diag);
end