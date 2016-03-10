function [CU] = curls2(X)

%
% curl operator
%
% X is a d-dimensional staggered grid s2
% CU is a d-dimensional staggered grid
%

d   = X.dim;

CU = staggered(d);

CU.M{1} = diff(X.M{3},[],2)*d(2)-diff(X.M{2},[],3)*d(3);
CU.M{2} = diff(X.M{1},[],3)*d(3)-diff(X.M{3},[],1)*d(1);
CU.M{3} = diff(X.M{2},[],1)*d(1)-diff(X.M{1},[],2)*d(2);

