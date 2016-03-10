function [Y] = curl_adjs2(dxS,dyS,dzS,z)

%
% adjoint curl operator
%
% z is a d-dimensional staggered grid
% Y is a d-dimensional staggered grid s2
%
%
d   = z.dim;

Y = staggereds2(d);

Y.M{1} = -(dyS(z.M{3})*d(2)-dzS(z.M{2})*d(3));
Y.M{2} = -(dzS(z.M{1})*d(3)-dxS(z.M{3})*d(1));
Y.M{3} = -(dxS(z.M{2})*d(1)-dyS(z.M{1})*d(2));