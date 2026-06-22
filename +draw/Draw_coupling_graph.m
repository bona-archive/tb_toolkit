% function [xvec, yvec] = Draw_coupling_graph(x,y,dist,linecolor,doPlot,color_distinguish_scatter)
%     if nargin < 5 || isempty(doPlot)
%         doPlot = true;
%     end
%     if nargin < 6 || isempty(color_distinguish_scatter)
%         do_color_distinguish_scatter = false;
%     else
%         do_color_distinguish_scatter = true;
%     end
%     ind1i=[];ind2i=[];
%     % NL = length(x);
%     x = x(:);
%     y = y(:);
%     for i=1:length(x)
%         ind1=1:length(x);
%         ind2=i*ones(1,length(x));
% 
%         dx=x-x(i);
%         dy=y-y(i);
%         dxy=(dx.^2+dy.^2).^(1/2); % calculate distance
% 
%         tempind=ind1(abs(dxy-dist)<1e-6);
%         ind1i=[ind1i ind1(tempind)];
%         ind2i=[ind2i ind2(tempind)];
% 
%     end
% 
%     xvec=reshape([x(ind1i) x(ind2i) x(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
%     yvec=reshape([y(ind1i) y(ind2i) y(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
%     %%%%%%%%% test draw %%%%%%%%%%%%%%%%%
%     if doPlot && do_color_distinguish_scatter
% 
%         % figure
%         line(xvec,yvec,'Color',linecolor);
%         hold on
%         scatter(x,y,[],color_distinguish_scatter,'filled')
%         axis equal
%     elseif doPlot && ~do_color_distinguish_scatter
%         line(xvec,yvec,'Color',linecolor);
%         hold on
%         scatter(x,y,'k','filled')
%         axis equal
%     else
%     end
% 
% 
% end

%%


function [] = Draw_coupling_graph(x,y,sublat,hops,doPlot,do_color_distinguish_scatter)
    if nargin < 5 || isempty(doPlot)
        doPlot = true;
    end
    if nargin < 6 || isempty(do_color_distinguish_scatter)
        do_color_distinguish_scatter = false;
    else
        do_color_distinguish_scatter = true;
    end
        % NL = length(x);
    x = x(:);
    y = y(:);
    sublat = sublat(:);

    edge_default = [
                    0.0, 0.0, 0.0;   % 검정
                    0.7, 0.35, 0.0;  % 갈색
                    0.4, 0.0, 0.7;   % 보라
                    0.0, 0.6, 0.5;   % 청록
                    ];
    node_default = [
                    1.0, 0.0, 0.0;   % 빨강
                    0.0, 0.0, 1.0;   % 파랑
                    0.0, 0.75, 0.0;  % 초록 (순수 [0,1,0]은 너무 형광)
                    1.0, 0.55, 0.0;  % 주황
                    0.6, 0.0, 0.8;   % 보라
                    0.0, 0.7, 0.7;   % 청록
                   ];

    t_vals = unique(hops(:,4));
    n_t = length(t_vals);
    edge_colors = edge_default(mod((1:n_t)-1, size(edge_default,1))+1, :);

    sub_vals = unique(sublat);
    n_sub = length(sub_vals);
    node_colors = node_default(mod((1:n_sub)-1, size(node_default,1))+1, :);
    
    for t_ind = 1:length(t_vals)
        t_val = t_vals(t_ind);
        ind1i=[];ind2i=[];
        for i=1:length(x)
            ind1=1:length(x);
            ind2=i*ones(1,length(x));
        
            dx=x-x(i);
            dy=y-y(i);
            dxy=(dx.^2+dy.^2).^(1/2); % calculate distance
            sub1 = sublat(i);
            for hop_ind = 1: size(hops,1)
                hop = hops(hop_ind,:);
                sub1_check = hop(1);
                sub2_check = hop(2);
                dist = hop(3);
                if sub1 == sub1_check && hop(4) == t_val
                    tempind=ind1(abs(dxy-dist)<1e-6 & sublat == sub2_check);
                    ind1i=[ind1i ind1(tempind)];
                    ind2i=[ind2i ind2(tempind)];
                end
            end
        
        end
        xvec=reshape([x(ind1i) x(ind2i) x(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
        yvec=reshape([y(ind1i) y(ind2i) y(ind2i)+NaN]',[3*length(x(ind2i)) 1]);


        %%%%%%%%% test draw %%%%%%%%%%%%%%%%%
        if doPlot
            hold on
            line(xvec,yvec,'Color',edge_colors(t_ind, :),'LineWidth',2);
        end

    end
    if doPlot
        if do_color_distinguish_scatter
            hold on
            scatter(x,y,[],node_colors(sublat, :),'filled')
        else
            hold on
            scatter(x,y,'k','filled')
        end
        axis equal
    end


end