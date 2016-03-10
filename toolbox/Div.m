function [div] = Div(phi)

div = diff(phi.M{1},[],1)*size(phi.M{1},1)+diff(phi.M{2},[],2)*size(phi.M{2},2)+diff(phi.M{3},[],3)*size(phi.M{3},3);

div = max(max(max(div)));