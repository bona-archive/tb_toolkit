function [super_hops] = primhop2superhop(prim_hops,superlat_2d_x,superlat_2d_y, latvec_pbc)


[super_basis, ~, prim_super_sublat_map] = obc.lat_2d2superbasis(superlat_2d_x,superlat_2d_y);

super_hops = [];

for hops_ind = 1: size(prim_hops,1)
    hop = prim_hops(hops_ind,:);

    % si = hop(1);
    % sj = hop(2);
    % dist_ij = hop(3);
    % t_ij = hop(4);
    % 
    % I_sites = find(prim_super_sublat_map(:,1) == si);
    % J_sites = find(prim_super_sublat_map(:,1) == sj);
    % for I_ind = 1:length(I_sites)
    %     I = I_sites(I_ind);
    %     for J_ind = 1:length(J_sites)
    %         J = J_sites(J_ind);
            for n1 = -1:1
                d = n1*latvec_pbc + super_basis(J,:) - super_basis(I,:);
                if abs(norm(d)-dist_ij) < 1e-4
                    super_hops = [super_hops;[I, J, dist_ij, t_ij]];
                else
                end
            end
    %     end
    % end
end
super_hops = unique(super_hops, 'rows');
end