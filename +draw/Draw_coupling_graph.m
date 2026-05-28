function [xvec, yvec] = Draw_coupling_graph(x,y,dist,color,doPlot)
    if nargin < 5 || isempty(doPlot)
        doPlot = true;
    end
    ind1i=[];ind2i=[];
    NL = length(x);
    x = x(:);
    y = y(:);
    for i=1:length(x)
        ind1=1:length(x);
        ind2=i*ones(1,length(x));
    
        dx=x-x(i);
        dy=y-y(i);
        dxy=(dx.^2+dy.^2).^(1/2); % calculate distance
    
        tempind=ind1(abs(dxy-dist)<1e-6);
        ind1i=[ind1i ind1(tempind)];
        ind2i=[ind2i ind2(tempind)];
    
    end
    
    xvec=reshape([x(ind1i) x(ind2i) x(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
    yvec=reshape([y(ind1i) y(ind2i) y(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
    %%%%%%%%% test draw %%%%%%%%%%%%%%%%%
    if doPlot
    
        % figure
        line(xvec,yvec,'Color',color);
        hold on
        scatter(x,y,'k','filled')
        axis equal
    end


end