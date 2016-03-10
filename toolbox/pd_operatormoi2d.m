function Y = pd_operatormoi2d(X, direction)

%
% pd_operator - linear operator for primal-dual scheme
%
%   Y = operator_pd(X, direction)
%
%   direction==+1 to compute K
%   direction==-1 to compute K^*
%
%   Copyright (c) 2012 Gabriel Peyre
%


if(direction == 1)
    
    dx = @(V)(size(X,1)-1)*diff(V,[],1);
    dy = @(V)(size(X,1)-1)*diff(V,[],2);

    intx = @(V)0.5*(V(2:end,:)+V(1:end-1,:));
    inty = @(V)0.5*(V(:,2:end)+V(:,1:end-1));
    
    % compute K
    Y(:,:,1) = intx(dy(X));
    Y(:,:,2) = inty(-dx(X));
    
else    
    
    dxS = @(V)size(X,1)*cat(1,-V(1,:), -diff(V,[],1), V(end,:));
    dyS = @(V)size(X,1)*cat(2,-V(:,1), -diff(V,[],2), V(:,end));
    intxS = @(V)0.5*cat(1,V(1,:),V(1:end-1,:)+V(2:end,:),V(end,:));
    intyS = @(V)0.5*cat(2,V(:,1),V(:,1:end-1)+V(:,2:end),V(:,end));

    % compute K^*
    Y = dyS(intxS(X(:,:,1))) - dxS(intyS(X(:,:,2)));
    
end

