function [hops] = make_hops(basis, latvecs, neighbor_order, t_vals, options)

arguments (Input)
    basis (:,2)
    latvecs (2,2)
    neighbor_order (1,:)
    t_vals (1,:)
    options.doPlot (1,1) logical = false
end
assert(size(neighbor_order,2)==size(t_vals,2), 'neighbor_order와 t_vals의 사이즈가 다릅니다.')


latvec_x = latvecs(1,:);
latvec_y = latvecs(2,:);

n_range=-1:1;
[x,y] = tb.lat_2d(length(n_range),length(n_range),latvec_x,latvec_y,basis);
D = tb.distmat(x(:), y(:));
unique_dist = unique(D(:));
[sorted_D, ~] = sort(unique_dist);

if options.doPlot
    figure
    for neighbor_order_ind = 1: length(neighbor_order)
        subplot(1,length(neighbor_order),neighbor_order_ind)
        hold on
        [~,~] = tb.Draw_coupling_graph(x,y,sorted_D(neighbor_order(neighbor_order_ind)+1),options.doPlot);
        axis equal
    end
else
end

[~,~,sub_SUBLAT] = ind2sub(size(x),1:length(D));

i=0;
hops = [];
for neighbor_order_ind = 1: length(neighbor_order)
    i=i+1;
    for a = 1: length(D)
        for b = 1: length(D)
    
            if abs(D(a,b)-sorted_D(neighbor_order(neighbor_order_ind)+1)) < 1e-6
                hops=[hops;sub_SUBLAT(a),sub_SUBLAT(b),sorted_D(neighbor_order(neighbor_order_ind)+1),t_vals(i)];
            else
    
            end
    
        end
    end
end
hops = unique(hops, 'rows');

end