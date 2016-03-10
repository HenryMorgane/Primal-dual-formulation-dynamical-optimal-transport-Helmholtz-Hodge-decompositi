function Y = pd_operator(dxS,dyS,dzS,X, direction)
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
    % compute K
    Y = interp(curls2(X));
else    
    % compute K^*
    Y = curl_adjs2(dxS,dyS,dzS,interp_adj(X));
end

